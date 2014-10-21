module Structure

  helper :mutated_isoforms_to_residue_list do |mutated_isoforms|
    residues = {}

    begin
      TSV.traverse Annotated.purge(mutated_isoforms), :type => :array do |mi|
        protein, _sep, change = mi.partition ":"
        if change.match(/([A-Z])(\d+)([A-Z])$/)
          next if $1 == $3
          position = $2
          residues[protein] ||=[]
          residues[protein] << position.to_i
        end
      end
    rescue Exception
      Log.exception $!
      raise $!
    end

    organism = mutated_isoforms.respond_to?(:organism)? mutated_isoforms.organism || Organism.default_code("Hsa") : Organism.default_code("Hsa")

    TSV.setup(residues, :key_field => "Ensembl Protein ID", :fields => ["Residues"], :type => :flat, :cast => :to_i, :namespace => organism)

    residues
  end

  helper :mutated_isoforms do |mutations,organism,watson|
    raise ParameterException, "No mutated_isoforms or genomic_mutations specified" if mutations.nil? 
    log :mutated_isoforms, "Finding mutated isoforms for genomic mutations"

    job = Sequence.job(:mutated_isoforms_fast, clean_name, :mutations => mutations, :organism => organism, :watson => watson)
    job.run true

    mis = Set.new
    TSV.traverse job do |m,_mis|
      mis.merge _mis 
    end
    mis = mis.to_a
    job.join

    FileUtils.mkdir_p files_dir
    if Open.remote? job.path
      Open.write(file(:mutated_isoforms_for_genomic_mutations), Open.read(job.path))
    else
      FileUtils.cp job.path, file(:mutated_isoforms_for_genomic_mutations)
    end

    mis
  end

  helper :residue_neighbours do |residues,organism|
    log :residue_neighbours, {:desc => "Find neighbouring residues"} do |bar|

      all_neighbours = TSV.setup({}, :key_field => "Isoform:residue", :fields => ["Ensembl Protein ID", "Residue", "PDB", "Neighbours"], :type => :list)

      TSV.traverse(residues, :bar => bar, :into => all_neighbours, :cpus => $cpus) do |protein, list|
        list = list.flatten.compact.uniq

        Structure.neighbours_i3d(protein, list, organism)
      end

      all_neighbours
    end
  end

  helper :isoform_to_mutation do |mi_annotations|
    log :mapping, "Mapping isoform annotations to genomic mutations"
    index = file(:mutated_isoforms_for_genomic_mutations).tsv :key_field => "Mutated Isoform", :merge => true, :type => :flat, :unnamed => true

    mi_parser = TSV::Parser.new mi_annotations

    mutation_annotation_dumper = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => ["Mutated Isoform"] + mi_parser.fields[1..-1], :type => :double, :namespace => mi_parser.namespace
    mutation_annotation_dumper.init
    TSV.traverse mi_parser, :into => mutation_annotation_dumper, :bar => "Isoform => mutation" do |mi, values|
      mi = mi.first if Array === mi
      mutations = index[mi]
      new_values = [[mi]] + values.collect{|v| [v * ";"] }

      mutations.each do |mutation|
        mutation_annotation_dumper.add mutation, new_values
      end
      nil
    end

    stream = mutation_annotation_dumper.stream

    log :collapsing, "Isoform to mutation stream collapsing"
    stream_collapsed = TSV.collapse_stream stream
    log :collapsed, "Isoform to mutation stream collapsed"

    Misc.save_stream file(:genomic_mutation_annotations), stream_collapsed.stream
  end

  helper :residue_to_isoform do |mis,residue_annotations|
    log :mapping, "Mapping residue annotations to isoforms"

    mi_annotation_dumper = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => residue_annotations.fields[1..-1], :type => :double, :namespace => residue_annotations.namespace 
    mi_annotation_dumper.init

    residue_annotations.unnamed = true
    TSV.traverse mis, :type => :array, :bar => {:desc => "Residue => Isoform", :max => mis.length},  :into => mi_annotation_dumper do |mi|
      next if mi.nil? or mi.empty?
      protein, change = mi.split(":")
      next if change.nil?
      
      if change.match(/([A-Z])(\d+)([A-Z])$/)
        begin
          next if $1 == $3
          position = $2.to_i
          raise "No Match" if residue_annotations[protein].nil?
          entries = Misc.zip_fields(residue_annotations[protein])
          entries.select!{|residue, *rest| residue.to_i == position}
          next if entries.empty?
          entries.each{|p| p.shift }

          [mi, Misc.zip_fields(entries.uniq)]
        rescue
          next
        end
      else
        next
      end
    end

    stream = mi_annotation_dumper.stream

    Misc.save_stream file(:mutated_isoform_annotations), stream
  end

  helper :residue_to_isoform_neighbours do |mis,residue_annotations,neighbours|
    log :mapping, "Mapping residue annotations to isoform neighbours"

    mi_annotation_dumper = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => residue_annotations.fields[1..-1], :type => :double, :namespace => residue_annotations.namespace 
    mi_annotation_dumper.init

    residue_annotations.unnamed = true
    TSV.traverse mis, :type => :array, :into => mi_annotation_dumper do |mi|
      next if mi.nil? or mi.empty?
      protein, change = mi.split(":")
      next if change.nil?
      if change.match(/([A-Z])(\d+)([A-Z])$/)
        next if $1 == $3
        mutated_position = $2
        isoforms_residue = [protein, mutated_position] * ":"
        next unless neighbours.include? isoforms_residue
        neighbour_residues = neighbours[isoforms_residue][-1].split ";"
        neighbour_residues.each do |position|
          begin
            position = position.to_i
            raise "No Match" if residue_annotations[protein].nil?
            entries = Misc.zip_fields(residue_annotations[protein])
            entries.select!{|residue, *rest| residue.to_i == position}
            next if entries.empty?
            entries.each{|p| p.shift }

            mi_annotation_dumper.add mi, Misc.zip_fields(entries.uniq)
          rescue
            next
          end
        end
      end
      nil
    end

    stream = mi_annotation_dumper.stream

    out, save = Misc.tee_stream stream

    Thread.new(Thread.current) do |parent|
      begin
        Misc.sensiblewrite(file(:mutated_isoform_annotations), save)
        save.join if save.respond_to? :join
      rescue
        save.abort if save.respond_to? :abort
        parent.raise $!
      end
    end

    out
  end
end
