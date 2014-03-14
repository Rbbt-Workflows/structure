module Structure
 
  desc <<-EOF
Given a set of proteins and resudies inside these proteins, finds the protein features that overlap, as annotated in UniProt.

The proteins and residues are specified as a TSV file, with Ensembl Protein IDs
as key, and residues as values. Types `flat` (separated by more tabs),
and `double` (separated by '|') are supported.
  EOF
  input :residues, :tsv, "Proteins and their affected residues", nil
  task :annotate_residues_UNIPROT => :tsv do |residues|
    raise ParameterException, "No residues provided" if residues.nil?
    tsv = TSV.setup({}, :key_field => "Ensembl Protein ID", :fields => ["Residue", "UniProt Features", "UniProt Feature locations", "UniProt Feature Descriptions"], :type => :double)

    iso2uni = Organism.protein_identifiers("Hsa").index :target => "UniProt/SwissProt Accession", :persist => true
    iso2sequence = Organism.protein_sequence("Hsa").tsv :type => :single, :persist => true

    missing = []
    residues.each do |isoform, list|
      list = Array === list ? list.flatten : [list]

      uniprot = iso2uni[isoform]
      if uniprot.nil?
        missing << isoform
        next
      end

      features = Structure.corrected_uniprot_features(uniprot, iso2sequence[isoform])
      overlapping = [[],[],[],[]]
      list.each do |position|
        position = position.to_i
        features.select{|info|
          case info[:type]
          when "VAR_SEQ", "CONFLICT", "CHAIN", "UNSURE"
            false
          when "DISULFID", "CROSSLNK", "VARIANT"
            info[:start] == position or info[:end] == position
          else
            info[:start] <= position and info[:end] >= position
          end
        }.each{|info|
          overlapping[0] << position
          overlapping[1] << info[:type]
          overlapping[2] << [info[:start], info[:end]] * ":"
          overlapping[3] << (info[:description] || "").strip.sub(/\.$/,'')
        }
      end

      tsv[isoform] = overlapping
    end

    if missing.any?
      Log.warn "Some isoforms failed to translate to UniProt: #{missing.length}"
      set_info(:missing, missing)
    end

    tsv
  end
  export_asynchronous :annotate_residues_UNIPROT

  desc <<-EOF
Given a set of proteins and resudies inside these proteins, finds the protein features that overlap, as annotated in Appris.

The proteins and residues are specified as a TSV file, with Ensembl Protein IDs
as key, and residues as values. Types `flat` (separated by more tabs),
and `double` (separated by '|') are supported.
  EOF
  input :residues, :tsv, "Proteins and their affected residues", nil
  task :annotate_residues_Appris => :tsv do |residues|
    tsv = TSV.setup({}, :key_field => "Ensembl Protein ID", :fields => ["Residue", "Appris Features", "Appris Feature locations", "Appris Feature Descriptions"], :type => :double)

    iso2sequence = Organism.protein_sequence("Hsa").tsv :type => :single, :persist => true

    residues.each do |isoform, list|
      list = Array === list ? list.flatten : [list]

      features = Structure.appris_features(isoform)

      overlapping = [[],[],[],[]]
      list.each do |position|
        position = position.to_i
        features.select{|info|
          info[:start] <= position and info[:end] >= position
        }.each{|info|
          overlapping[0] << position
          overlapping[1] << info[:type]
          overlapping[2] << [info[:start], info[:end]] * ":"
          overlapping[3] << (info[:description] || "").strip.sub(/\.$/,'')
        }
      end

      tsv[isoform] = overlapping
    end

    tsv
  end
  export_asynchronous :annotate_residues_Appris

  desc <<-EOF
Given a set of proteins and resudies inside these proteins, finds the mutations registered in COSMIC that affect those residues, and provide
some annotations from the samples that contained them

The proteins and residues are specified as a TSV file, with Ensembl Protein IDs
as key, and residues as values. Types `flat` (separated by more tabs),
and `double` (separated by '|') are supported.

The result is the proteins along with the overlapping features and some information about them
  EOF
  input :residues, :tsv, "Proteins and their affected residues", nil
  task :annotate_variants_COSMIC => :tsv do |residues|

    cosmic_residue_mutations = Structure.COSMIC_residues

    isoform_matched_variants = {}
    residues.through do |protein,list|
      list = Array === list ? list.flatten : [list]

      list.each do |position|
        matching_mutations = cosmic_residue_mutations[[protein, position]*":"]
        next if matching_mutations.nil? or matching_mutations.empty?
        isoform_matched_variants[protein] ||= []
        isoform_matched_variants[protein].concat matching_mutations
      end
    end

    cosmic_mutation_annotations = Structure.COSMIC_mutation_annotations

    res = TSV.setup({}, :key_field => "Ensembl Protein ID", :fields => ["Residue", "Genomic Mutation"].concat(cosmic_mutation_annotations.fields), :type => :double)

    isoform_matched_variants.each do |protein, mutations|
      values = []
      mutations.each do |mutation|
        begin
          annotations = cosmic_mutation_annotations[mutation]
          raise "No annotations for #{ mutation } in #{cosmic_mutation_annotations.persistence_path}" if annotations.nil?
          aa_mutation = (annotations.first || []).compact.first
          raise "No AA mutation" if aa_mutation.nil?
          if aa_mutation.match(/(\d+)/)
            residue = $1
          else
            raise "AA mutation does not have a position: #{ aa_mutation }" if aa_mutation.nil?
          end
          values << [residue,mutation].concat(annotations || [])
        rescue
          Log.exception $!
          next
        end
      end

      res[protein] = Misc.zip_fields(values)
    end

    res
  end
  export_asynchronous :annotate_variants_COSMIC


  desc <<-EOF
Given a set of proteins and resudies inside these proteins, finds the known variants that overlap, as annotated in UniProt.

The proteins and residues are specified as a TSV file, with Ensembl Protein IDs
as key, and residues as values. Types `flat` (separated by more tabs),
and `double` (separated by '|') are supported.

The result is the proteins along with the overlapping features and some information about them
  EOF
  input :residues, :tsv, "Proteins and their affected residues", nil
  task :annotate_variants_UNIPROT => :tsv do |residues|

    uniprot_residue_mutations = Structure.UniProt_residues

    isoform_matched_variants = {}
    residues.each do |protein, list|
      list = Array === list ? list.flatten : [list]


      list.each do |position|
        key = [protein, position]*":"

        matching_mutations = uniprot_residue_mutations[key]
        next if matching_mutations.nil? or matching_mutations.empty?
        isoform_matched_variants[protein] ||= []
        isoform_matched_variants[protein].concat matching_mutations
      end
    end

    uniprot_mutation_annotations = Structure.UniProt_mutation_annotations

    res = TSV.setup({}, :key_field => "Ensembl Protein ID", :fields => ["Residue"] + uniprot_mutation_annotations.fields, :type => :double)

    isoform_matched_variants.each do |protein, mutations|
      values = []
      mutations.each do |mutation|
        begin
          annotations = uniprot_mutation_annotations[mutation]
          raise "No annotations for #{ mutation } in #{uniprot_mutation_annotations.persistence_path}" if annotations.nil?
          residue = annotations.first.scan(/\d+/)
          values << [residue].concat(annotations || [])
        rescue
          Log.exception $!
        end
      end

      res[protein] = Misc.zip_fields(values)
    end

    res
  end
  export_asynchronous :annotate_variants_UNIPROT


  input :residues, :tsv, "Proteins and their affected residues", nil
  task :residue_interfaces => :tsv do |residues|

    neighbour_residues = {}

    residues.with_monitor :desc => "Processing residues" do
    residues.through do |protein, list|
      list = list.flatten.compact.sort
      begin
        neighbours = Structure.interface_neighbours_i3d(protein.dup, list)
        next if neighbours.nil? or neighbours.empty?
        neighbour_residues.merge!(neighbours)
      rescue
        Log.exception $!
        next
      end
    end
    end

    TSV.setup(neighbour_residues, :key_field => "Isoform:residue:partner", :fields => ["Partner", "PDB", "Partner residues"], :type => :double)
  end
  export_asynchronous :residue_interfaces

end
