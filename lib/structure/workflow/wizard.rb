
module Structure
  input :mutations, :array, "Mutations (e.g. 18:6237978:G, ENSP00000382976:L257R, L3MBTL4:L257R)"
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  task :wizard => :tsv do |mutations,organism|
    databases = Structure::ANNOTATORS.keys 

    raise ParameterException, "No mutations given" if mutations.nil? 
    raise ParameterException, "This wizard is limited to a thousand variants. For larger jobs use standard functions please" if mutations.length > 1000

    log :identifiying
    mutations = mutations.collect{|m| m.sub(/(.*) (.*)/,'\1:\2')}
    type = case mutations.first
           when nil
             raise ParameterException, "No mutations given"
           when /.{1,2}:(\d+):[ACTG+-]+/
             :genomic
           when /ENS.*:\w+\d+\w+/
             :protein
           when /.*:\w+\d+\w+/
             index = Organism.identifiers(organism).index :target => "Ensembl Gene ID", :order => true, :persist => true
             ensg2enst = Organism.transcripts(organism).tsv :key_field => "Ensembl Gene ID", :fields => ["Ensembl Transcript ID"], :persist => true, :merge => true, :type => :flat
             enst2ensp = Organism.transcripts(organism).index :target => "Ensembl Protein ID", :fields => ["Ensembl Transcript ID"], :persist => true, :merge => true
             enst2name = Organism.transcript_name(organism).index :target => "Ensembl Transcript Name", :fields => ["Ensembl Transcript ID"], :persist => true, :merge => true
             uni2ensp = Organism.uniprot2ensembl(organism).index :persist => true, :target => "Ensembl Protein ID"
             ensp2uni = Organism.ensembl2uniprot(organism).index :persist => true, :target => "UniProt/SwissProt Accession"
             mutations = mutations.collect do |m|
                orig_gene, _sep, change = m.partition ":"
                gene = index[orig_gene]
                gene_transcripts = ensg2enst[gene].sort_by{|t| enst2name[t].split("-").last.to_i}

                gene_isoforms = enst2ensp.values_at *gene_transcripts
                perfect_isoforms = gene_isoforms.select{|p| 
                  _uni = ensp2uni[p]
                  _p = uni2ensp[_uni]
                  p == _p
                }

                principal_transcripts = (Appris::PRINCIPAL_TRANSCRIPTS & gene_transcripts)
                principal_isoforms = enst2ensp.values_at *principal_transcripts
                uni_pricipal_isoforms = principal_isoforms.select{|p| ensp2uni[p]}

                perfect_principal_isoforms = uni_pricipal_isoforms & perfect_isoforms

                if perfect_principal_isoforms.any?
                  protein = perfect_principal_isoforms.first
                elsif uni_pricipal_isoforms.any?
                  protein = uni_pricipal_isoforms
                elsif perfect_isoforms.any?
                  protein = perfect_isoforms.first
                elsif principal_isoforms.any?
                  protein = principal_isoforms.first
                else
                  protein = gene_isoforms.first
                end

                [protein, change] * ":"
             end
             :protein
           end

    case type
    when :genomic
      all_annotations = Sequence.job(:mutated_isoforms_fast, "Wizard", :mutations => mutations, :organism => organism, :principal => true).run.to_double
      databases.each do |database|
        log database
        annotations = Structure.job(:annotate, "Wizard", :mutations => mutations, :organism => organism, :database => database, :principal => true).run
        Open.write(file(database), annotations.to_s)
        all_annotations = all_annotations.attach(annotations)
      end

      log :interfaces
      interfaces = Structure.job(:interfaces, "Wizard", :mutations => mutations, :organism => organism).run
      interfaces.rename_field "Ensembl Protein ID", "Partner Ensembl Protein ID"
      Open.write(file('interfaces'), interfaces.to_s)
      all_annotations = all_annotations.attach(interfaces)

      all_annotations_n = Sequence.job(:mutated_isoforms_fast, "Wizard", :mutations => mutations, :organism => organism, :principal => true).run.to_double
      databases.each do |database|
        log database + ' neighbours'
        annotations = Structure.job(:annotate_neighbours, "Wizard", :mutations => mutations, :organism => organism, :database => database, :principal => true).run
        Open.write(file(database + ' neighbours'), annotations.to_s)
        annotations.rename_field "Residue", database + " residue"
        all_annotations_n = all_annotations_n.attach(annotations)
      end

    when :protein

      all_annotations = TSV.setup(mutations, :key_field => "Mutated Isoform", :fields => [], :type => :double, :namespace => organism)
      databases.each do |database|
        log database
        annotations = Structure.job(:annotate_mi, "Wizard", :mutated_isoforms => mutations, :organism => organism, :database => database).run
        Open.write(file(database), annotations.to_s)
        all_annotations = all_annotations.attach(annotations)
      end

      log :interfaces
      interfaces = Structure.job(:mi_interfaces, "Wizard", :mutated_isoforms => mutations, :organism => organism).run
      interfaces.rename_field "Ensembl Protein ID", "Partner Ensembl Protein ID"
      Open.write(file('interfaces'), interfaces.to_s)
      all_annotations = all_annotations.attach(interfaces)

      all_annotations_n = TSV.setup(mutations, :key_field => "Mutated Isoform", :fields => [], :type => :double, :namespace => organism)
      databases.each do |database|
        log database + ' neighbours'
        annotations = Structure.job(:annotate_mi_neighbours, "Wizard", :mutated_isoforms => mutations, :organism => organism, :database => database).run
        Open.write(file(database + ' neighbours'), annotations.to_s)
        annotations.rename_field "Residue", database + " residue"
        all_annotations_n = all_annotations_n.attach(annotations)
      end
    end

    all_annotations_n.fields = all_annotations_n.fields.collect{|f| "Neighbour " + f }
    all_annotations.attach(all_annotations_n, :fields => all_annotations_n.fields.reject{|f| f =~ /Mutated Isoform/})

    all_annotations
  end
  export_asynchronous :wizard

  def self.score_for(field, value, all_values = nil)
    value = value.collect{|v| v.split(";")}.flatten
    score = case field
            when "Appris Features"
              if value.include? "firestar"
                2
              else
                1 
              end
            when "UniProt Features"
              relevant = %w(DISULFID DNA_BIND METAL INTRAMEM CROSSLNK MUTAGEN)
              value = value.zip(all_values["UniProt Feature Descriptions"].collect{|v| v.split(";")}.flatten).reject{|v,d| v == "MUTAGEN" and d =~ /no effect/i}.collect{|v,d| v}
              sum = 0
              sum += 1 if (value & relevant).any?
              sum
            when "Sample ID"
              case 
              when value.length > 10
                3
              when value.length > 5
                2
              when value.length > 1
                1 
              else
                0
              end
            when "UniProt Variant ID"
              case 
              when value.length > 0
                1
              else
                0
              end
            when "Type of Variant"
              if value.include?("Disease")
                2
              elsif value.include?("Unclassified")
                1
              else
                0
              end
            when "Partner Ensembl Protein ID"
              2
            else
              0
            end
    score
  end

  def self.score_mi(values)
    score = 0
    values.zip(values.fields).each do |value, field|
      next if value.empty?
      if field =~ /Neighbour/
        score = score.to_f + (score_for(field.sub('Neighbour ',''), value, values).to_f / 2)
      else
        score = score + score_for(field, value, values)
      end
    end
    score
  end

  dep :wizard
  task :scores => :tsv do 
    wizard = step(:wizard).load

    wizard.add_field "Score" do |mi, values|
      Structure.score_mi(values)
    end

    wizard.reorder "Mutated Isoform", "Score"
  end
  export_asynchronous :scores
end
