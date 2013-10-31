require 'rbbt/sources/uniprot'
require 'rbbt/sources/COSMIC'
require 'rbbt/workflow'

Workflow.require_workflow "Appris"

module Structure
  local_persist_setup

  def self.COSMIC_mutations_over_protein(protein, organism = "Hsa/jan2013")
    mutated_isoforms = local_persist_tsv(COSMIC.mutations, "Mutated Isoforms") do |data|
      mutations = local_persist("COSMIC Mutation") do
        pos = COSMIC.mutations.open do |file|
          TSV.parse_header(file).fields.index("Genomic Mutation")
        end
        CMD.cmd("grep -v '^#'|cut -f#{pos+2}", :in => COSMIC.mutations.open).read.split("\n").reject{|l| l.empty?}.uniq.sort
      end

      job = Sequence.job(:mutated_isoforms_for_genomic_mutations, "Cosmic", :mutations => mutations, :organism => organism, :watson => true)
      job.run(true).path.tsv(:persist => true, :persist_data => data)

      data
    end

    mutated_isoforms.with_unnamed do
      mutated_isoforms.with_monitor do
        mutated_isoforms.select("Mutated Isoform" => /^#{protein}:[A-Z]\d+[A-Z*]$/)
      end
    end
  end
 
  input :gene, :string, "Ensembl Gene ID"
  input :organism, :string, "Organism code", "Hsa/jan2013"
  task :gene_principal_isoform => :string do |gene, organism|
    gene = Gene.setup(gene.dup, "Ensembl Gene ID", organism)
    pis = gene.principal_isoforms
    if pis.nil? or pis.empty?
      log :isoform, "Using longest isoform"
      gene.proteins.sort_by{|protein| protein.sequence.length}.last
    else
      log :isoform, "Using principal isoform"
      pis.first.protein
    end
  end

  dep :gene_principal_isoform
  input :organism, :string, "Organism code", "Hsa/jan2013"
  task :gene_mutations => :tsv do |organism|
    protein = step(:gene_principal_isoform).load
    mutations = Structure.COSMIC_mutations_over_protein(protein, organism)
    changes = TSV.setup({}, :key_field => "Genomic Mutation", :fields => ["Change"], :type => :single)
    log :changes, "Determining mutation changes in isoform"
    mutations.with_unnamed do
      mutations.through do |mutation, isoforms|
        changes[mutation] = isoforms.select{|mi| mi =~  /^#{protein}:[A-Z]\d+[A-Z*]/}.collect{|mi| mi.split(":").last}.reject{|c| c[0] == c[-1]}.first
      end
    end
    changes
  end

  dep :gene_mutations
  input :distance, :float, "Distance in Armstrongs", 4
  input :organism, :string, "Organism code", "Hsa/jan2013"
  task :gene_mutation_close_residues => :tsv do |distance, organism|
    gene_mutations = step(:gene_mutations).path.tsv :type => :double

    isoform = step(:gene_principal_isoform).load
    Protein.setup(isoform, "Ensembl Protein ID", organism)
    sequence = isoform.sequence

    pdb_scores = {}
    isoform_pdbs = isoform.pdbs || {}
    isoform_pdbs.each{|pdb, info|
      if info[:region].nil?
        pdb_scores[pdb] = 0
      else
        pdb_scores[pdb] = (info[:region].end - info[:region].begin) * (info[:method] == "X-ray" ? 1.1 : 1)
      end
    }

    pdbs = isoform_pdbs.keys

    pdb_aligments = {}
    pdb_close_contacts = {}
    log :pdb, "Aligning pdbs to isoform sequence and determining close contacts"
    pdbs.each do |pdb|
      pdb_aligments[pdb] = Structure.job(:alignment_map, isoform, :sequence => sequence, :pdb => pdb).run(true).path.tsv :type => :list
      pdb_close_contacts[pdb] = PDBHelper.pdb_close_residues(distance, pdb)
    end
    pdb_inverse_aligments = {}
    pdb_aligments.each do |pdb,tsv|
      pdb_inverse_aligments[pdb] = tsv.reorder(tsv.fields.first)
    end

    log :close_residues, "Finding residues close to mutation sites"
    gene_mutations.add_field "Close residues" do |mutation, values|
      changes = values["Change"]
      position = changes.collect{|c| c.match(/(\d+)/)[1].to_i}.first.to_s
      best_pdb = pdb_aligments.select{|pdb,alignment|
        alignment[position]
      }.sort_by{|pdb,alignment|
        pdb_scores[pdb]
      }.collect{|pdb,a| pdb}.last

      if best_pdb.nil?
        []
      else
        pdb_aligments[best_pdb][position].collect do |pdb_pos|
          close = pdb_close_contacts[best_pdb][pdb_pos]
          next if close.nil?
          close.collect{|p| 
            pdb_inverse_aligments[best_pdb][p]
          }.compact
        end.compact.flatten
      end
    end

    gene_mutations
  end

  dep :gene_mutation_close_residues
  input :organism, :string, "Organism code", "Hsa/jan2013"
  task :gene_mutation_uniprot_info => :tsv do |organism|
    require 'rbbt/sources/uniprot'

    gene_mutation_close_residues = step(:gene_mutation_close_residues).load
    isoform = step(:gene_principal_isoform).load
    Protein.setup(isoform, "Ensembl Protein ID", organism)
    
    uniprot = isoform.uniprot

    raise "No UniProt for the principal isoform" if uniprot.nil?

    log :alignment_map, "Aligning UniProt sequence to isoform sequence"
    uniprot_alignment, ensembl_alignment = SmithWaterman.align(UniProt.sequence(uniprot), isoform.sequence)
    alignment_map = Structure.alignment_map(uniprot_alignment, ensembl_alignment)

    aligned_features = {}
    UniProt.features(uniprot).each do |info|
      type, start, eend, description = info.values_at :type, :start, :end, :description
      start = alignment_map[start]
      eend = alignment_map[eend]

      next if start.nil? or eend.nil?

      (start.to_i..eend.to_i).to_a.collect do |pos|
        aligned_features[pos] = {:type => type, :description => description}
      end
    end

    log :features, "Finding matching features"
    gene_mutation_close_residues.add_field "UniProt features" do |mutation, values|
      values["Change"].collect{|change|
        position = change.match(/(\d+)/)[1].to_i
        aligned_features[position.to_i]
      }.compact.uniq.collect{|info| [info[:type], info[:description]] * ": "}
    end

    log :features, "Finding close features"
    gene_mutation_close_residues.add_field "Neighbour UniProt features" do |mutation, values|
      values["Close residues"].collect{|pos|
        aligned_features[pos.to_i]
      }.compact.uniq.collect{|info| [info[:type], info[:description]] * ": "}
    end

    gene_mutation_close_residues
  end

  export_asynchronous :gene_mutation_uniprot_info

  dep :gene_mutation_close_residues
  input :organism, :string, "Organism code", "Hsa/jan2013"
  task :gene_mutation_appris_info => :tsv do |organism|
    require 'rbbt/sources/uniprot'

    gene_mutation_close_residues = step(:gene_mutation_close_residues).load
    isoform = step(:gene_principal_isoform).load
    Protein.setup(isoform, "Ensembl Protein ID", organism)
    
    aligned_features = {}
    isoform.appris_residues.each do |info|
      type, ranges = info
      ranges.each do |v|
        start, eend = v.values_at "start", "end"
        description = [start, eend] * ".."

        (start.to_i..eend.to_i).to_a.collect do |pos|
          aligned_features[pos] = {:type => type, :description => description}
        end
      end
    end

    log :features, "Finding matching features"
    gene_mutation_close_residues.add_field "UniProt features" do |mutation, values|
      values["Change"].collect{|change|
        position = change.match(/(\d+)/)[1].to_i
        aligned_features[position.to_i]
      }.compact.uniq.collect{|info| [info[:type], info[:description]] * ": "}
    end

    log :features, "Finding close features"
    gene_mutation_close_residues.add_field "Neighbour UniProt features" do |mutation, values|
      values["Close residues"].collect{|pos|
        aligned_features[pos.to_i]
      }.compact.uniq.collect{|info| [info[:type], info[:description]] * ": "}
    end

    gene_mutation_close_residues
  end

  export_asynchronous :gene_mutation_appris_info

  dep :gene_mutations
  task :gene_mutation_phenotype => :tsv do
    gene_mutations = step(:gene_mutations).load
    cosmic_info = COSMIC.mutations.tsv(:persist => true, :key_field => "Genomic Mutation", :fields => ["Primary site", "Site subtype", "Primary histology", "Histology subtype"], :type => :double, :merge => true)

    gene_mutations.attach cosmic_info, :fields => cosmic_info.fields
  end

  export_asynchronous :gene_mutation_phenotype


end
