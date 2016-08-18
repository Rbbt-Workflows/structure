
module Structure
  input :mutations, :array, "Mutations (e.g. 18:6237978:G, ENSP00000382976:L257R, L3MBTL4:L257R)"
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  input :watson, :boolean, "Mutations all reported on the watson (forward) strand as opposed to the gene strand", true
  task :wizard => :tsv do |mutations,organism,watson|
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
                next if gene.nil? or ensg2enst[gene].nil? 
                gene_transcripts = ensg2enst[gene].sort_by{|t| enst2name[t].split("-").last.to_i}

                gene_isoforms = enst2ensp.values_at(*gene_transcripts).compact.reject{|i| i.empty?}
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

                next if protein.nil? or protein.empty?

                [protein, change] * ":"
             end.compact
             :protein
           end

    case type
    when :genomic
      all_annotations = Sequence.job(:mutated_isoforms_fast, "Wizard", :mutations => mutations, :organism => organism, :principal => true, :watson => watson).run.to_double
      databases.each do |database|
        log database
        annotations = Structure.job(:annotate, "Wizard", :mutations => mutations, :organism => organism, :database => database, :principal => true, :watson => watson).run
        Open.write(file(database), annotations.to_s)
        all_annotations = all_annotations.attach(annotations)
      end

      log :interfaces
      interfaces = Structure.job(:interfaces, "Wizard", :mutations => mutations, :organism => organism).run
      interfaces.rename_field "Ensembl Protein ID", "Partner Ensembl Protein ID"
      Open.write(file('interfaces'), interfaces.to_s)
      all_annotations = all_annotations.attach(interfaces)

      all_annotations_n = Sequence.job(:mutated_isoforms_fast, "Wizard", :mutations => mutations, :organism => organism, :principal => true, :watson => watson).run.to_double
      databases.each do |database|
        log database + ' neighbours'
        annotations = Structure.job(:annotate_neighbours, "Wizard", :mutations => mutations, :organism => organism, :database => database, :principal => true, :watson => watson).run
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

    all_annotations.namespace = organism

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
    wizard = step(:wizard)
    wizard_res = wizard.load

    wizard_res.add_field "Score" do |mi, values|
      Structure.score_mi(values)
    end

    wizard_res.reorder "Mutated Isoform", "Score"
  end
  
  dep :scores
  task :score_summary => :tsv do 
    Step.wait_for_jobs dependencies
    scores = step(:scores)
    wizard = scores.step(:wizard)
    wizard_res = wizard.load
    scores_res = scores.load


    require 'rbbt/sources/clinvar'
    Workflow.require_workflow "DbNSFP"
    Workflow.require_workflow "Pandrugs"


    clinvar = ClinVar.mi_summary.tsv :fields => ["ClinicalSignificance"], :type => :single, :persist => true
    interfaces = wizard.file('interfaces').tsv
    dbNSFP = DbNSFP.job(:annotate, nil, :mutations => wizard_res.keys).run
    cosmic = wizard.file('COSMIC').tsv
    cosmic_neighbours = wizard.file('COSMIC neighbours').tsv
    appris = wizard.file('Appris').tsv
    appris_neighbours = wizard.file('Appris neighbours').tsv
    uniprot = wizard.file('UniProt').tsv
    uniprot_neighbours = wizard.file('UniProt neighbours').tsv
    predictors = %w(SIFT Polyphen2_HDIV Polyphen2_HVAR MutationTaster MutationAssessor FATHMM LRT VEST3 CADD )
    thresholds = %w( <0.05 >0.957,0.453 >0.909,0.447 >0.5 >3.5,1.9 <-1.5 - >0.8 >3.5   )

    fields = %w(Score CV #CS FL MR MUT DT PPI DP)
    tsv = TSV.setup({}, :key_field => "Mutated Isoform", :fields => fields, :type => :list, :namespace => wizard_res.namespace)

    scores_res.each do |mutation, score|
      values = []
      ensp, _sep, change = mutation.partition(":") 
      next unless ensp =~ /^ENSP/

      damage_count = 0
      total_preds = 0
      dvalues = dbNSFP[mutation]
      if dvalues
        predictors.each_with_index do |predictor,i|
          next if predictor == "LRT"
          raw, dscore, converted, rankscore, raw_rankscore, converted_rankscore, p = nil
          threshold = thresholds[i]
          raw = dvalues[predictor + '_raw'] if dvalues.fields.include? predictor + '_raw'
          dscore = dvalues[predictor + '_score'] if dvalues.fields.include? predictor + '_score'
          dscore = nil if String === dscore and dscore.empty?
          dscore = raw if dscore.nil?
          converted = dvalues[predictor + '_converted_score'] if dvalues.fields.include? predictor + '_converted_score'
          rankscore = dvalues[predictor + '_rankscore'] if dvalues.fields.include? predictor + '_rankscore'
          raw_rankscore = dvalues[predictor + '_raw_rankscore'] if dvalues.fields.include? predictor + '_raw_rankscore'
          converted_rankscore = dvalues[predictor + '_converted_rankscore'] if dvalues.fields.include? predictor + '_converted_rankscore'

          if score and threshold != '-'
            p = case threshold
              when /^<(.*)/
                ths = $1.split(",")
                ths.inject(0){|acc,e| acc += 1 if dscore.to_f < e.to_f; acc}.to_f/ths.length
              when /^>(.*)/
                ths = $1.split(",")
                ths.inject(0){|acc,e| acc += 1 if dscore.to_f > e.to_f; acc}.to_f/ths.length
              else
                nil
              end

            damage_count += 1 if p > 0.5
            total_preds +=1

          end
        end
      end

      values << score.flatten.first

      if clinvar.include? mutation and clinvar[mutation] == "Pathogenic"
        values << "Yes"
      else
        values << "No"
      end

      if cosmic.include? mutation
        count = cosmic[mutation]["Sample ID"].length
      else
        count = 0
      end

      if cosmic_neighbours.include? mutation
        ncount = cosmic_neighbours[mutation]["Sample ID"].collect{|l| l.split(";")}.flatten.uniq.length
      else
        ncount = 0
      end

      values << "#{count} (#{ncount})"


      if appris.include? mutation
        count = Misc.zip_fields(appris[mutation]).select{|type,lig| type =~ /firestar/}.length
      else
        count = 0
      end

      if appris_neighbours.include? mutation
        ncount = Misc.zip_fields(appris_neighbours[mutation]).select{|res,type,lig| type =~ /firestar/}.length
      else
        ncount = 0
      end
      values << "#{count == 0 ? "No" : "Yes"} (#{ncount})"



      if uniprot.include? mutation
        count = Misc.zip_fields(uniprot[mutation]).select{|feat,loc,desc| feat =~ /MOD_RES/}.length
      else
        count = 0
      end

      if uniprot_neighbours.include? mutation
        ncount = Misc.zip_fields(uniprot_neighbours[mutation]).select{|res,feat,loc,desc| feat =~ /MOD_RES/}.length
      else
        ncount = 0
      end

      values << "#{count == 0 ? "No" : "Yes"} (#{ncount})"

      if uniprot.include? mutation
        count = Misc.zip_fields(uniprot[mutation]).select{|feat,loc,desc| feat =~ /MUTAGEN/}.length
      else
        count = 0
      end

      if uniprot_neighbours.include? mutation
        ncount = Misc.zip_fields(uniprot_neighbours[mutation]).select{|res,feat,loc,desc| feat =~ /MUTAGEN/}.length
      else
        ncount = 0
      end

      values << "#{count == 0 ? "No" : "Yes"} (#{ncount})"

      if Pandrugs.knowledge_base.subset('gene_drugs', :source => [mutation.protein.gene], :target => :all).filter(:target_marker => 'target').filter(:status => "Approved").length > 0
        values << "Yes"
      else
        values << "No"
      end

      values << interfaces.include?(mutation) ? "Yes" : "No"

      if dbNSFP.include? mutation
        values << "#{damage_count} of 8"
      else
        values << "NA"
      end
    
      tsv[mutation] = values
    end

    tsv
  end

  export_asynchronous :scores, :score_summary
end
