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

    organism = mutated_isoforms.respond_to?(:organism)? mutated_isoforms.organism || "Hsa" : "Hsa"

    TSV.setup(residues, :key_field => "Ensembl Protein ID", :fields => ["Residues"], :type => :flat, :cast => :to_i, :namespace => organism)

    residues
  end

  helper :mutated_isoforms do |mutations,organism,watson|
    raise ParameterException, "No mutated_isoforms or genomic_mutations specified" if mutations.nil? 
    log :mutated_isoforms, "Finding mutated isoforms for genomic mutations: #{ watson.inspect }"

    job = Sequence.job(:mutated_isoforms, clean_name, :mutations => mutations, :organism => organism, :watson => watson)
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
    log :residue_neighbours, "Find neighbouring residues"

    all_neighbours = TSV.setup({}, :key_field => "Isoform:residue", :fields => ["Ensembl Protein ID", "Residue", "PDB", "Neighbours"], :type => :list)

    TSV.traverse(residues, :into => all_neighbours, :cpus => $cpus) do |protein, list|
      list = list.flatten.compact.uniq
      Structure.neighbours_i3d(protein, list, organism)
    end

    all_neighbours
  end

  helper :residue_to_isoform do |mis,residue_annotations|
    log :mapping, "Mapping residue annotations to isoforms"

    mi_annotation_dumper = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => residue_annotations.fields[1..-1], :type => :double, :namespace => residue_annotations.namespace 
    mi_annotation_dumper.init

    TSV.traverse mis, :type => :array, :cpus => 2, :into => mi_annotation_dumper do |mi|
      next if mi.nil? or mi.empty?
      protein, change = mi.split(":")
      next if change.nil?
      
      if change.match(/([A-Z])(\d+)([A-Z])$/)
        begin
          next if $1 == $3
          position = $2.to_i
          raise "No Match" if residue_annotations[protein].nil?
          entries = residue_annotations[protein].zip_fields
          entries.select!{|residue, *rest| residue.to_i == position}
          next if entries.empty?
          entries.each{|p| p.shift }

          [mi, Misc.zip_fields(entries.uniq)]
        rescue
          Log.exception $!
          next
        end
      else
        next
      end
    end

    stream = mi_annotation_dumper.stream

    out, save = Misc.tee_stream stream

    Thread.new do
      Misc.sensiblewrite(file(:mutated_isoform_annotations), save)
    end

    out
  end

  helper :isoform_to_mutation do |mi_annotations|
    log :mapping, "Mapping isoform annotations to genomic mutations"
    index = file(:mutated_isoforms_for_genomic_mutations).tsv :key_field => "Mutated Isoform", :merge => true, :type => :flat

    mi_parser = TSV::Parser.new mi_annotations

    mutation_annotation_dumper = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => ["Mutated Isoform"] + mi_parser.fields[1..-1], :type => :double, :namespace => mi_parser.namespace
    mutation_annotation_dumper.init
    TSV.traverse mi_parser, :into => mutation_annotation_dumper do |mi, values|
      mi = mi.first if Array === mi
      mutations = index[mi]
      mutations.each do |mutation|
        new_values = [[mi]] + values.collect{|v| [v * ";"] }
        mutation_annotation_dumper.add mutation, new_values
      end
      nil
    end

    stream = mutation_annotation_dumper.stream

    stream_collapsed = TSV.collapse_stream stream

    out, save = Misc.tee_stream stream_collapsed.stream

    Thread.new do
      Misc.sensiblewrite(file(:genomic_mutation_annotations), save)
    end

    out
  end
end
