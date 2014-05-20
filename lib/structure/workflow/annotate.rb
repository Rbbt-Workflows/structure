require 'structure/annotator'
module Structure

  CORRECTED_FEATURES = Rbbt.var.cache.Structure.corrected_features.find
  Open.repository_dirs << CORRECTED_FEATURES unless Open.repository_dirs.include? CORRECTED_FEATURES

  NEIGHBOURS = Rbbt.var.cache.Structure.neighbours.find
  Open.repository_dirs << NEIGHBOURS unless Open.repository_dirs.include? NEIGHBOURS

  INTERFACE_NEIGHBOURS = Rbbt.var.cache.Structure.interface_neighbours.find
  Open.repository_dirs << INTERFACE_NEIGHBOURS unless Open.repository_dirs.include? INTERFACE_NEIGHBOURS

  ANNOTATORS = IndiferentHash.setup({})

  ANNOTATORS["COSMIC"] = Annotator.new "Genomic Mutation", 'Sample ID', 'Primary site', 'Site subtype', 'Primary histology', 'Histology subtype',
    'Genome-wide screen', 'Mutation somatic status', 'Sample source', 'Tumour origin' do |isoform, residue,organism|

    @cosmic_residue_mutations ||= Structure.COSMIC_residues
    @cosmic_mutation_annotations ||= Structure.COSMIC_mutation_annotations
    isoform_residue = isoform + ":" << residue.to_s
    mutations = @cosmic_residue_mutations[isoform_residue]
    next if mutations.nil? 
    mutations.uniq!
    next if mutations.empty?
    tmp = {}
    mutations.each do |mutation|
      annot = @cosmic_mutation_annotations[mutation]
      Misc.zip_fields(annot).each do |a|
        sample, *rest = a
        next if tmp.include? sample
        tmp[sample] = rest
      end
    end
    annot = tmp.collect{|p| p.flatten}
    annot = Misc.zip_fields(annot)
    annot.unshift mutations
    annot
  end

  ANNOTATORS["Appris"] = Annotator.new "Appris Features", "Appris Feature locations", "Appris Feature Descriptions" do |isoform, residue,organism|
    features = Structure.appris_features(isoform)

    overlapping = [[],[],[]]
    features.select{|info|
      info[:start] <= residue and info[:end] >= residue
    }.each{|info|
      overlapping[0] << info[:type]
      overlapping[1] << [info[:start], info[:end]] * ":"
      overlapping[2] << (info[:description] || "").strip.sub(/\.$/,'')
    }

    next if overlapping.first.empty?
    overlapping
  end

  ANNOTATORS["InterPro"] = Annotator.new "InterPro ID", "Domain range" do |isoform,residue,organism|
    @iso2uni ||= {}
    @iso2sequence ||= {}
    iso2uni = @iso2uni[organism] ||= Organism.protein_identifiers(organism).index(:target => "UniProt/SwissProt Accession", :persist => true)
    iso2sequence = @iso2sequence[organism] ||= Organism.protein_sequence(organism).tsv(:type => :single, :persist => true)

    uniprot = iso2uni[isoform]
    next if uniprot.nil?
    sequence = iso2sequence[isoform]
    next if sequence.nil?

    features =  Persist.persist("Corrected InterPro features", :yaml, :dir => CORRECTED_FEATURES, :other => {:uniprot => uniprot, :sequence => sequence}) do 
      Structure.corrected_interpro_features(uniprot, sequence)
    end
    next if features.empty?

    overlapping = [[],[]]
    features.select{|info|
      info[:start] <= residue and info[:end] >= residue
    }.each{|info|
      overlapping[0] << info[:code]
      overlapping[1] << [info[:start], info[:end]] * ":"
    }

    next if overlapping.first.empty?
    overlapping
  end

  ANNOTATORS["UniProt"] = Annotator.new "UniProt Features", "UniProt Feature locations", "UniProt Feature Descriptions" do |isoform,residue,organism|
    @iso2uni ||= {}
    @iso2sequence ||= {}
    iso2uni = @iso2uni[organism] ||= Organism.protein_identifiers(organism).index(:target => "UniProt/SwissProt Accession", :persist => true)
    iso2sequence = @iso2sequence[organism] ||= Organism.protein_sequence(organism).tsv(:type => :single, :persist => true)

    uniprot = iso2uni[isoform]
    next if uniprot.nil?
    sequence = iso2sequence[isoform]
    next if sequence.nil?

    features =  Persist.persist("Corrected UniProt features", :yaml, :dir => CORRECTED_FEATURES, :other => {:uniprot => uniprot, :sequence => sequence}) do 
      Structure.corrected_uniprot_features(uniprot, sequence)
    end
    next if features.empty?

    overlapping = [[],[],[]]

    features.select{|info|
      case info[:type]
      when "VAR_SEQ", "CONFLICT", "CHAIN", "UNSURE"
        false
      when "DISULFID", "CROSSLNK", "VARIANT"
        info[:start] == residue or info[:end] == residue
      else
        info[:start].to_i <= residue and info[:end].to_i >= residue
      end
    }.each{|info|
      description = (info[:description] || "").strip.sub(/\.$/,'')
      #description = "-" if description.nil? or description.empty?
      overlapping[0] << info[:type]
      overlapping[1] << [info[:start], info[:end]] * ":"
      overlapping[2] << description
    }

    next if overlapping.first.empty?

    overlapping
  end

  ANNOTATORS["variants"] = Annotator.new "UniProt Variant ID", "SNP ID",  "Type of Variant", "Disease" do |isoform,residue,organism|
    @iso2uni ||= {}
    @iso2sequence ||= {}
    @annotations ||= Structure.UniProt_mutation_annotations
    iso2uni = @iso2uni[organism] ||= Organism.protein_identifiers(organism).index(:target => "UniProt/SwissProt Accession", :persist => true)
    iso2sequence = @iso2sequence[organism] ||= Organism.protein_sequence(organism).tsv(:type => :single, :persist => true)

    uniprot = iso2uni[isoform]
    next if uniprot.nil?
    sequence = iso2sequence[isoform]
    next if sequence.nil?

    features =  Persist.persist("Corrected UniProt features", :yaml, :dir => CORRECTED_FEATURES, :other => {:uniprot => uniprot, :sequence => sequence}) do 
      Structure.corrected_uniprot_features(uniprot, sequence)
    end
    next if features.empty?

    overlapping = [[],[],[],[]]

    features.select{|info|
      info[:type] == "VARIANT" and info[:start] == residue
    }.each{|info|
      if info[:description].match(/(VAR_\d+)/)
        id = $1
        next unless @annotations.include? id
        overlapping[0] << id
        annots = @annotations[id]
        overlapping[1] << annots[2]
        overlapping[2] << annots[1]
        overlapping[3] << annots[3]
      end
    }

    next if overlapping.first.empty?

    overlapping
  end

  #{{{ NEIGHBOURS
  
  input :mutated_isoforms, :array, "Mutated Isoform", nil, :stream => true
  input :organism, :organism, "Organism code", "Hsa"
  task :mi_neighbours => :tsv do |mis,organism|
    annotations = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => ["Residue", "PDB", "Neighbours"], :type => :double
    annotations.init
    TSV.traverse mis, :cpus => $cpus, :bar => "Mutated Isoform neighbours", :into => annotations, :type => :array do |mi|
      case
      when mi =~ /^(.*):([A-Z])(\d+)([A-Z])$/
        next if $2 == $4
        isoform = $1
        residue = $3.to_i
      when mi =~ /^(.*):(\d+)$/
        isoform = $1
        residue = $2.to_i
      else
        next
      end

      n =  Persist.persist("Neighbours", :yaml, :dir => NEIGHBOURS, :other => {:isoform => isoform, :residue => residue, :organism => organism}) do 
        Structure.neighbours_i3d(isoform, [residue], organism)
      end

      next if n.empty?
      pdbs = []
      ns = []
      n.each do |r,v|
        _p,_pos,pdb,n = v
        next if n.empty?
        pdbs << pdb
        ns << n.split(";")
      end
      next if ns.empty?
      [mi, [residue, pdbs, ns]]
    end
  end

  input :mutated_isoforms, :array, "Mutated Isoform", nil, :stream => true
  input :organism, :string, "Organism code", "Hsa"
  task :mi_interfaces => :tsv do |mis,organism|
    mi_annotations = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => ["Residue", "Partner Ensembl Protein ID", "PDB", "Partner Residues"], :type => :double, :namespace => organism
    mi_annotations.init
    TSV.traverse mis, :cpus => $cpus, :bar => "Mutated Isoform neighbours", :into => mi_annotations, :type => :array do |mi|
      case
      when mi =~ /^(.*):([A-Z])(\d+)([A-Z])$/
        next if $2 == $4
        isoform = $1
        residue = $3.to_i
      when mi =~ /^(.*):(\d+)$/
        isoform = $1
        residue = $2.to_i
      else
        next
      end

      n = Persist.persist("Interface neighbours", :yaml, :dir => INTERFACE_NEIGHBOURS, :other => {:isoform => isoform, :residue => residue, :organism => organism}) do 
        Structure.interface_neighbours_i3d(isoform.dup, [residue], organism)
      end

      next if n.nil? or n.empty?

      all_annots = []
      n.each do |r,v|
        _pos,part,pdb,n = v
        next if part.nil? or part.empty?
        all_annots << [residue, part, pdb, n]
      end

      [mi, Misc.zip_fields(all_annots)]
    end
  end

  input :mutated_isoforms, :array, "Mutated Isoform", nil, :stream => true
  input :organism, :string, "Organism code", "Hsa"
  input :database, :select, "Database of annotations", "UniProt", :select_options => ANNOTATORS.keys
  task :annotate_mi => :tsv do |mis,organism,database|
    annotator = ANNOTATORS[database]
    raise ParameterException, "Database not identified: #{ database }" if annotator.nil?
    annotator.organism = organism

    mi_annotations = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => ["Residue"].concat(annotator.fields), :type => :double, :namespace => organism
    mi_annotations.init
    TSV.traverse mis, :cpus => $cpus, :bar => "Annot. #{ database }", :into => mi_annotations, :type => :array do |mi|
      case
      when mi =~ /^(.*):([A-Z])(\d+)([A-Z])$/
        next if $2 == $4
        isoform = $1
        residue = $3.to_i
      when mi =~ /^(.*):(\d+)$/
        isoform = $1
        residue = $2.to_i
      else
        next
      end

      annotations = annotator.annotate isoform, residue, organism
      next if annotations.nil?

      annotations.unshift residue

      [mi, annotations]
    end
  end

  dep :mi_neighbours
  input :database, :select, "Database of annotations", "UniProt", :select_options => ANNOTATORS.keys
  task :annotate_mi_neighbours => :tsv do |database|
    neigh = step(:mi_neighbours).grace
    organism = neigh.info[:inputs][:organism]

    annotator = ANNOTATORS[database]
    raise ParameterException, "Database not identified: #{ database }" if annotator.nil?
    annotator.organism = organism

    mi_annotations = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => ["Residue"].concat(annotator.fields), :type => :double, :namespace => organism
    mi_annotations.init
    TSV.traverse neigh, :into => mi_annotations  do |mi,v|
      mi = mi.first if Array === mi
      case
      when mi =~ /^(.*):([A-Z])(\d+)([A-Z])$/
        next if $2 == $4
        isoform = $1
        residue = $3.to_i
      when mi =~ /^(.*):(\d+)$/
        isoform = $1
        residue = $2.to_i
      else
        next
      end

      res, pdb, residues = v

      all_annots = []
      residues.each do |residue|
        residue = residue.to_i

        annotations = annotator.annotate isoform, residue, organism
        next if annotations.nil?

        annotations.unshift residue
        all_annots << annotations.collect{|v| Array === v ? v*";" : v }
      end

      next if all_annots.empty?
      [mi, Misc.zip_fields(all_annots.uniq)]
    end
  end

  #{{{ MUTATION ANNOTATIONS

  dep Sequence, :mutated_isoforms_fast
  input :database, :select, "Database of annotations", "UniProt", :select_options => ANNOTATORS.keys
  task :annotate => :tsv do |database|
    mutated_isoforms = step(:mutated_isoforms_fast)
    mutated_isoforms.grace
    organism = mutated_isoforms.info[:inputs][:organism]

    annotator = ANNOTATORS[database]
    raise ParameterException, "Database not identified: #{ database }" if annotator.nil?
    annotator.organism = organism

    mutation_annotations = TSV::Dumper.new :key_field => "Genomic Mutations", :fields => ["Mutated Isoform", "Residue"].concat(annotator.fields), :type => :double, :namespace => organism
    mutation_annotations.init

    TSV.traverse mutated_isoforms, :cpus => $cpus, :bar => "Annot. #{ database }", :into => mutation_annotations do |mutation,mis|
      next if mis.nil? or mis.empty?
      all_annots = []

      mis.each do |mi|
        case
        when mi =~ /^(.*):([A-Z])(\d+)([A-Z])$/
          next if $2 == $4
          isoform = $1
          residue = $3.to_i
        when mi =~ /^(.*):(\d+)$/
          isoform = $1
          residue = $2.to_i
        else
          next
        end

        annotations = annotator.annotate isoform, residue, organism
        next if annotations.nil?
        annotations = annotations.collect{|v| v * ";" }

        annotations.unshift residue
        annotations.unshift mi

        all_annots << annotations
      end

      next if all_annots.empty?
      [mutation, Misc.zip_fields(all_annots)]
    end
  end

  dep Sequence, :mutated_isoforms_fast
  input :database, :select, "Database of annotations", "UniProt", :select_options => ANNOTATORS.keys
  task :annotate_neighbours => :tsv do |database|
    mutated_isoforms = step(:mutated_isoforms_fast)
    mutated_isoforms.grace
    organism = mutated_isoforms.info[:inputs][:organism]

    annotator = ANNOTATORS[database]
    raise ParameterException, "Database not identified: #{ database }" if annotator.nil?
    annotator.organism = organism

    mutation_annotations = TSV::Dumper.new :key_field => "Genomic Mutations", :fields => ["Mutated Isoform", "Residue"].concat(annotator.fields), :type => :double, :namespace => organism
    mutation_annotations.init

    TSV.traverse mutated_isoforms, :cpus => $cpus, :bar => "Annot. #{ database }", :into => mutation_annotations do |mutation,mis|
      next if mis.nil? or mis.empty?
      all_annots = []

      used_mi = nil
      mis.each do |mi|
        case
        when mi =~ /^(.*):([A-Z])(\d+)([A-Z])$/
          next if $2 == $4
          isoform = $1
          residue = $3.to_i
        when mi =~ /^(.*):(\d+)$/
          isoform = $1
          residue = $2.to_i
        else
          next
        end

        n =  Persist.persist("Neighbours", :yaml, :dir => NEIGHBOURS, :other => {:isoform => isoform, :residue => residue, :organism => organism}) do 
          Structure.neighbours_i3d(isoform, [residue], organism)
        end
        next if n.empty?

        pdbs = []
        ns = []
        n.each do |r,v|
          _p,_pos,pdb,n = v
          next if n.empty?
          pdb
          n.split(";").each do |nresidue|

            annotations = annotator.annotate isoform, nresidue.to_i, organism
            next if annotations.nil?
            annotations = annotations.collect{|v| v * ";" }

            annotations.unshift nresidue

            all_annots << annotations
          end
        end
        used_mi = mi
        break
      end

      next if all_annots.empty?
      [mutation, [used_mi] + Misc.zip_fields(all_annots)]
    end
  end

  dep Sequence, :mutated_isoforms_fast
  task :interfaces => :tsv do |mis,organism|
    mutated_isoforms = step(:mutated_isoforms_fast)
    mutated_isoforms.grace
    organism = mutated_isoforms.info[:inputs][:organism]

    annotations = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => ["Partner Ensembl Protein ID", "PDB", "Partner Residues"], :type => :double, :namespace => organism
    annotations.init
    TSV.traverse mutated_isoforms, :cpus => $cpus, :bar => "Mutated Isoform neighbours", :into => annotations do |mutation,mis|
      next if mis.nil? or mis.empty?
      all_annots = []
      mis.each do |mi|
        case
        when mi =~ /^(.*):([A-Z])(\d+)([A-Z])$/
          next if $2 == $4
          isoform = $1
          residue = $3.to_i
        when mi =~ /^(.*):(\d+)$/
          isoform = $1
          residue = $2.to_i
        else
          next
        end

        n = Persist.persist("Interface neighbours", :yaml, :dir => INTERFACE_NEIGHBOURS, :other => {:isoform => isoform, :residue => residue, :organism => organism}) do 
          Structure.interface_neighbours_i3d(isoform.dup, [residue], organism)
        end

        next if n.nil? or n.empty?

        n.each do |r,v|
          _pos,part,pdb,n = v
          next if part.nil? or part.empty?
          all_annots << [part, pdb, n]
        end
      end

      next if all_annots.empty?
      [mutation, Misc.zip_fields(all_annots.uniq)]
    end
  end

  export_asynchronous :mi_neighbours, :mi_interfaces, :annotate_mi, :annotate_mi_neighbours, :interfaces, :annotate, :annotate_neighbours
end
