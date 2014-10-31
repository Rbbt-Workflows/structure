
module Structure
  input :mutations, :array, "Mutations (e.g. 18:6237978:G, ENSP00000382976:L257R, L3MBTL4:L257R)"
  input :organism, :string, "Organism code", Organism.default_code("Hsa")
  task :wizard => :tsv do |mutations,organism|
    databases = Structure::ANNOTATORS.keys 

    raise ParameterException, "This wizard is limited to a thousand variants. For larger jobs use standard functions please" if mutations.length > 1000

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
             gene2isoform = Organism.transcripts(organism).tsv :key_field => "Ensembl Gene ID", :fields => ["Ensembl Protein ID"], :type => :flat, :persist => true, :merge => true
             Log.tsv gene2isoform
             mutations = mutations.collect do |m|
                orig_gene, _sep, change = m.partition ":"
                gene = index[orig_gene]
                gene_proteins = gene2isoform[gene]
                protein = ( Appris::PRINCIPAL_ISOFORMS & gene_proteins).first
                [protein, change] * ":"
             end
             :protein
           end

    case type
    when :genomic
      all_annotations = Sequence.job(:mutated_isoforms_fast, "Wizard", :mutations => mutations, :organism => organism, :principal => true).run.to_double
      databases.each do |database|
        annotations = Structure.job(:annotate, "Wizard", :mutations => mutations, :organism => organism, :database => database, :principal => true).run
        Open.write(file(database), annotations.to_s)
        all_annotations = all_annotations.attach(annotations)
      end
      interfaces = Structure.job(:interfaces, "Wizard", :mutations => mutations, :organism => organism).run
      interfaces.rename_field "Ensembl Protein ID", "Partner Ensembl Protein ID"
      Open.write(file('interfaces'), interfaces.to_s)
      all_annotations = all_annotations.attach(interfaces)

      all_annotations_n = Sequence.job(:mutated_isoforms_fast, "Wizard", :mutations => mutations, :organism => organism, :principal => true).run.to_double
      databases.each do |database|
        annotations = Structure.job(:annotate_neighbours, "Wizard", :mutations => mutations, :organism => organism, :database => database, :principal => true).run
        Open.write(file(database + ' neighbours'), annotations.to_s)
        all_annotations_n = all_annotations_n.attach(annotations)
      end

    when :protein
      all_annotations = TSV.setup(mutations, :key_field => "Mutated Isoform", :fields => [], :type => :double, :namespace => organism)
      databases.each do |database|
        annotations = Structure.job(:annotate_mi, "Wizard", :mutated_isoforms => mutations, :organism => organism, :database => database).run
        Open.write(file(database), annotations.to_s)
        all_annotations = all_annotations.attach(annotations)
      end
      interfaces = Structure.job(:mi_interfaces, "Wizard", :mutated_isoforms => mutations, :organism => organism).run
      interfaces.rename_field "Ensembl Protein ID", "Partner Ensembl Protein ID"
      Open.write(file('interfaces'), interfaces.to_s)
      all_annotations = all_annotations.attach(interfaces)

      all_annotations_n = TSV.setup(mutations, :key_field => "Mutated Isoform", :fields => [], :type => :double, :namespace => organism)
      databases.each do |database|
        annotations = Structure.job(:annotate_mi_neighbours, "Wizard", :mutated_isoforms => mutations, :organism => organism, :database => database).run
        Open.write(file(database + ' neighbours'), annotations.to_s)
        all_annotations_n = all_annotations_n.attach(annotations)
      end
    end

    all_annotations_n.fields = all_annotations_n.fields.collect{|f| "Neighbour " + f }
    all_annotations.attach(all_annotations_n, :fields => all_annotations_n.fields.reject{|f| f =~ /Mutated Isoform/})

    all_annotations
  end

  export_asynchronous :wizard
end
