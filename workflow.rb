require 'rbbt'
require 'rbbt/util/misc'
require 'rbbt/workflow'

require 'rbbt/sources/organism'
require 'rbbt/sources/uniprot'
require 'rbbt/sources/InterPro'

require 'structure/alignment'
require 'structure/pdb_alignment'

require 'structure/pdb_helper'
require 'structure/neighbours'
require 'structure/interactome_3d'

require 'structure/uniprot'
require 'structure/appris'
require 'structure/COSMIC'
require 'structure/interpro'

Workflow.require_workflow 'Genomics'
Workflow.require_workflow 'Translation'
Workflow.require_workflow 'PdbTools'
Workflow.require_workflow "Sequence"

require 'rbbt/entity'
require 'rbbt/entity/mutated_isoform'
require 'structure/entity'

require 'rbbt/util/simpleopt'
$cpus ||= SOPT.get("--cpus*")[:cpus]
$cpus = $cpus.to_i if $cpus

Log.info "Loading Structure with #{ $cpus.inspect }" unless $cpus.nil?

module Structure
  extend Workflow

  input :genomic_mutations, :array, "Genomic Mutations"
  input :mutated_isoforms, :array, "Protein Mutations"
  input :database, :select, "Database of annotations", "UniProt", :select_options => ["UniProt", "COSMIC", "Appris", "InterPro", "variants"]
  input :organism, :string, "Organism code", "Hsa"
  input :watson, :boolean, "Mutations all reported on the watson (forward) strand as opposed to the gene strand", true
  task :annotated_variants => :tsv do |mutations, mis, database, organism,watson|
    raise ParameterException, "Database missing" unless database
    if mis.nil?
      raise ParameterException, "No genomic_mutations or mutated_isoforms" if mutations.nil?
      mis = mutated_isoforms mutations, organism, watson
    end

    residues = mutated_isoforms_to_residue_list(mis)

    log :annotating_residues, "Annotating residues with #{ database }"
    residue_annotations = case database
                          when "UniProt"
                            Structure.job(:annotate_residues_UniProt, clean_name, :residues => residues)
                          when "COSMIC"
                            Structure.job(:annotate_residues_COSMIC, clean_name, :residues => residues)
                          when "Appris"
                            Structure.job(:annotate_residues_Appris, clean_name, :residues => residues)
                          when "InterPro"
                            Structure.job(:annotate_residues_InterPro, clean_name, :residues => residues)
                          when "variants"
                            Structure.job(:annotate_residues_UniProt_variants, clean_name, :residues => residues)
                          else
                            raise ParameterException, "Unknown database: #{ Misc.fingerprint database }" 
                          end

    residue_annotations.clean if residue_annotations.error? or residue_annotations.aborted?

    residue_annotations.run false
    residue_annotations.grace
    stream = residue_annotations.get_stream || residue_annotations.path.open

    Misc.save_stream(file(:residue_annotations), stream)

    unless residue_annotations.done?
      residue_annotations.fork
      residue_annotations.join
    end

    mi_annotations = residue_to_isoform mis, residue_annotations.path

    if file(:mutated_isoforms_for_genomic_mutations).exists?
      isoform_to_mutation mi_annotations
    else
      mi_annotations
    end
  end
  export_asynchronous :annotated_variants


  input :genomic_mutations, :array, "Genomic Mutations"
  input :mutated_isoforms, :array, "Protein Mutations"
  input :database, :select, "Database of annotations", "UniProt", :select_options => ["UniProt", "COSMIC", "Appris", "InterPro", "variants"]
  input :organism, :string, "Organism code", "Hsa"
  input :watson, :boolean, "Mutations all reported on the watson (forward) strand as opposed to the gene strand", true
  task :annotated_variant_neighbours => :tsv do |mutations, mis, database, organism,watson|
    if mis.nil?
      mis = mutated_isoforms mutations, organism,watson
    end

    residues = mutated_isoforms_to_residue_list(mis)

    neighbours = self.residue_neighbours residues, organism

    Open.write(file(:neighbours), neighbours.to_s)

    neighbour_residues = {}

    log :neighbour_residues, "Sorting residues"
    neighbours.through do |iso,values|
      protein, _pos = iso.split ":"
      neighbouring_positions = values.last.split ";"

      neighbour_residues[protein] ||= []
      neighbour_residues[protein].concat neighbouring_positions
      neighbour_residues[protein].uniq!
    end
    residues.annotate neighbour_residues


    log :annotating_residues, "Annotating residues with #{ database }"
    residue_annotations = case database
                  when "InterPro"
                    Structure.job(:annotate_residues_InterPro, clean_name, :residues => neighbour_residues).fork
                  when "UniProt"
                    Structure.job(:annotate_residues_UniProt, clean_name, :residues => neighbour_residues).fork
                  when "variants"
                    Structure.job(:annotate_residues_UniProt_variants, clean_name, :residues => neighbour_residues).fork
                  when "COSMIC"
                    Structure.job(:annotate_residues_COSMIC, clean_name, :residues => neighbour_residues).fork
                  when "Appris"
                    Structure.job(:annotate_residues_Appris, clean_name, :residues => neighbour_residues).fork
                  end

    residue_annotations.join

    FileUtils.mkdir_p self.files_dir
    FileUtils.cp residue_annotations.path, file(:residue_annotations)

    Open.write(file(:residue_annotations), residue_annotations.to_s)

    mi_annotations = residue_to_isoform_neighbours mis, TSV.open(residue_annotations), neighbours

    if file(:mutated_isoforms_for_genomic_mutations).exists?
      isoform_to_mutation mi_annotations
    else
      mi_annotations
    end
  end
  export_asynchronous :annotated_variant_neighbours

  input :genomic_mutations, :array, "Genomic Mutations"
  input :mutated_isoforms, :array, "Protein Mutations"
  input :organism, :string, "Organism code", "Hsa"
  input :watson, :boolean, "Mutations all reported on the watson (forward) strand as opposed to the gene strand", true
  task :_variant_interfaces => :tsv do |mutations, mis, organism,watson|
    if mis.nil?
      mis = mutated_isoforms mutations, organism, watson
    end

    residues = mutated_isoforms_to_residue_list(mis)

    job = Structure.job(:residue_interfaces, clean_name, :residues => residues)
    residue_annotations = job.run

    Open.write(file(:residue_annotations), residue_annotations.to_s)

    mi_annotations = {}

    mis.each do |mi|
      protein, change = mi.split(":")
      if change.match(/([A-Z])(\d+)([A-Z])$/)
        begin
          next if $1 == $3
          position = $2.to_i
          raise "No Match" if residue_annotations[protein].nil?
          entries = residue_annotations[protein].zip_fields
          entries.select!{|residue, *rest| residue.to_i == position}
          next if entries.empty?
          entries.each{|p| p.shift }

          if mi_annotations[mi].nil?
            fixed_entries = Misc.zip_fields(entries.uniq)
            mi_annotations[mi] = fixed_entries
          else
            entries += Misc.zip_fields(mi_annotations[mi])
            fixed_entries = Misc.zip_fields(entries.uniq)
            mi_annotations[mi] = fixed_entries
          end
        rescue
          next
        end
      end
    end

    fields = residue_annotations.nil? ? [] : residue_annotations.fields[1..-1]

    TSV.setup(mi_annotations, :key_field => "Mutated Isoform", :fields => fields, :type => :double, :namespace => organism)

    mi_annotations

    Open.write(file(:mutated_isoform_annotations), mi_annotations.to_s)

    if file(:mutated_isoforms_for_genomic_mutations).exists?
      index = file(:mutated_isoforms_for_genomic_mutations).tsv :key_field => "Mutated Isoform", :merge => true, :type => :flat
      mutation_annotations = {}
      mi_annotations.each do |mi, values|
        mutations = index[mi]
        mutations.each do |mutation|
          new_values = [mi] + values.collect{|v| v * ";" }
          if mutation_annotations[mutation].nil?
            mutation_annotations[mutation] = new_values.collect{|v| [v] }
          else
            e = Misc.zip_fields(mutation_annotations[mutation])
            n = e << new_values
            mutation_annotations[mutation] = Misc.zip_fields(n)
          end
        end
      end

      TSV.setup(mutation_annotations, :key_field => "Genomic Mutation", :fields => ["Mutated Isoform"] + fields, :type => :double, :namespace => organism)
      Open.write(file(:genomic_mutation_annotations), mutation_annotations.to_s)

      mutation_annotations
    else
      mi_annotations
    end
  end
  export_asynchronous :variant_interfaces
end

require 'structure/workflow/alignments'
require 'structure/workflow/residues'
require 'structure/workflow/helpers'
require 'structure/workflow/annotate'
