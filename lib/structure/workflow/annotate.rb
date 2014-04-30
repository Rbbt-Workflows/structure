module Structure

  CORRECTED_FEATURES = Rbbt.var.cache.Structure.corrected_features.find
  Open.repository_dirs << CORRECTED_FEATURES unless Open.repository_dirs.include? CORRECTED_FEATURES

  NEIGHBOURS = Rbbt.var.cache.Structure.neighbours.find
  Open.repository_dirs << NEIGHBOURS unless Open.repository_dirs.include? NEIGHBOURS

  INTERFACE_NEIGHBOURS = Rbbt.var.cache.Structure.interface_neighbours.find
  Open.repository_dirs << INTERFACE_NEIGHBOURS unless Open.repository_dirs.include? INTERFACE_NEIGHBOURS

  ANNOTATORS = IndiferentHash.setup({})

  ANNOTATORS["COSMIC"] = Annotator.new "Genomic Mutation" do |isoform, residue|
    isoform_residue = isoform + ":" << residue.to_s
    mutations = cosmic_residue_mutations[isoform_residue]
    next if mutations.nil? 
    mutations.uniq!
    next if mutations.empty?
    annot = cosmic_mutation_annotations.values_at(*mutations)
    Misc.zip_fields(annot.compact).collect{|v| v * "|" }
    annot.unshift mutations * "|"
    annot
  end

  ANNOTATORS["Appris"] = Annotator.new "Appris Features", "Appris Feature locations", "Appris Feature Descriptions" do |isoform, residue|
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

  ANNOTATORS["InterPro"] = Annotator.new "InterPro ID", "Range" do |isoform,residue|
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

  ANNOTATORS["UniProt"] = Annotator.new "UniProt Features", "UniProt Feature locations", "UniProt Feature Descriptions" do |isoform,residue|
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
        info[:start] <= residue and info[:end] >= residue
      end
    }.each{|info|
      overlapping[0] << info[:code]
      overlapping[1] << [info[:start], info[:end]] * ":"
      overlapping[2] << (info[:description] || "").strip.sub(/\.$/,'')
    }

    next if overlapping.first.empty?

    overlapping
  end

  ANNOTATORS["variants"] = Annotator.new "UniProt Features", "UniProt Feature locations", "UniProt Feature Descriptions", "UniProt Variant ID", "SNP ID",  "Type of Variant", "Disease" do |isoform,residue|
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
      info[:type] == "VARIANT"
    }.each{|info|
      if info[:description].match(/(VAR_\d+)/)
        id = $1
        next unless annotations.include? id
        overlapping[0] << id
        annots = annotations[id]
        overlapping[1] << annots[2]
        overlapping[2] << annots[1]
        overlapping[3] << annots[3]
      end
    }

    next if overlapping.first.empty?

    overlapping
  end
  fields = ["UniProt Variant ID", "SNP ID",  "Type of Variant", "Disease"]
end

dep Sequence, :mutated_isoforms_fast
input :database, :select, "Database of annotations", "UniProt", :select_options => ["UniProt", "COSMIC", "Appris", "InterPro", "variants"]
task :annotate => :tsv do |database|

  mutated_isoforms = step(:mutated_isoforms_fast)
  mutated_isoforms.grace
  organism = mutated_isoforms.info[:inputs][:organism]

  annotator = ANNOTATORS[database]

  case database
  when "COSMIC"
    cosmic_residue_mutations = Structure.COSMIC_residues
    cosmic_mutation_annotations = Structure.COSMIC_mutation_annotations
    fields = ["Genomic Mutation"] + cosmic_mutation_annotations.fields
  when "variants"
    iso2uni = Organism.protein_identifiers(organism).index :target => "UniProt/SwissProt Accession", :persist => true
    iso2sequence = Organism.protein_sequence(organism).tsv :type => :single, :persist => true
    annotations = Structure.UniProt_mutation_annotations
    fields = ["UniProt Variant ID", "SNP ID",  "Type of Variant", "Disease"]
  when "UniProt", "variants"
    iso2uni = Organism.protein_identifiers(organism).index :target => "UniProt/SwissProt Accession", :persist => true
    iso2sequence = Organism.protein_sequence(organism).tsv :type => :single, :persist => true
    fields = ["UniProt Features", "UniProt Feature locations", "UniProt Feature Descriptions"]
  when "InterPro"
    iso2uni = Organism.protein_identifiers(organism).index :target => "UniProt/SwissProt Accession", :persist => true
    iso2sequence = Organism.protein_sequence(organism).tsv :type => :single, :persist => true
    fields = ["InterPro ID", "Range"]
  when "Appris"
    fields = ["Appris Features", "Appris Feature locations", "Appris Feature Descriptions"]
  end

  mutation_annotations = TSV::Dumper.new :key_field => "Genomic Mutations", :fields => ["Mutated Isoform", "Residue"].concat(fields), :type => :double, :namespace => organism
  mutation_annotations.init
  TSV.traverse mutated_isoforms, :cpus => $cpus, :bar => "Annot. #{ database }", :into => mutation_annotations do |mutation,mis|
    next if mis.nil? or mis.empty?
    all_mutation_annotations = []
    mis.each do |mi|
      next unless mi =~ /(.*):([A-Z])(\d+)([A-Z])$/
        next if $2 == $4
      isoform = $1
      residue = $3.to_i

      annotations = case database
                    when "COSMIC"
                      isoform_residue = isoform + ":" << residue.to_s
                      mutations = cosmic_residue_mutations[isoform_residue]
                        next if mutations.nil? 
                        mutations.uniq!
                        next if mutations.empty?
                        annot = cosmic_mutation_annotations.values_at(*mutations)
                        Misc.zip_fields(annot.compact).collect{|v| v * "|" }
                        annot.unshift mutations * "|"

                      when "Appris"
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

                      when "InterPro"
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

                      when "UniProt", "variants"
                        uniprot = iso2uni[isoform]
                        next if uniprot.nil?
                        sequence = iso2sequence[isoform]
                        next if sequence.nil?

                        features =  Persist.persist("Corrected UniProt features", :yaml, :dir => CORRECTED_FEATURES, :other => {:uniprot => uniprot, :sequence => sequence}) do 
                          Structure.corrected_uniprot_features(uniprot, sequence)
                        end
                        next if features.empty?

                        if database == "variants"
                          overlapping = [[],[],[],[]]
                        else
                          overlapping = [[],[],[]]
                        end
                        features.select{|info|
                          if database == "variants"
                            info[:type] == "VARIANT"
                          else
                            case info[:type]
                            when "VAR_SEQ", "CONFLICT", "CHAIN", "UNSURE"
                              false
                            when "DISULFID", "CROSSLNK", "VARIANT"
                              info[:start] == residue or info[:end] == residue
                            else
                              info[:start] <= residue and info[:end] >= residue
                            end
                          end
                        }.each{|info|
                          if database == "variants"
                            if info[:description].match(/(VAR_\d+)/)
                              id = $1
                              next unless annotations.include? id
                              overlapping[0] << id
                              annots = annotations[id]
                              overlapping[1] << annots[2]
                              overlapping[2] << annots[1]
                              overlapping[3] << annots[3]
                            end
                          else
                            overlapping[0] << info[:code]
                            overlapping[1] << [info[:start], info[:end]] * ":"
                            overlapping[2] << (info[:description] || "").strip.sub(/\.$/,'')
                          end
                        }

                        next if overlapping.first.empty?

                        overlapping
                      end

        annotations.unshift residue
        annotations.unshift mi
        all_mutation_annotations << annotations.collect{|v| Array === v ? v * ";" : v }
      end

      next if all_mutation_annotations.empty?
      [mutation, Misc.zip_fields(all_mutation_annotations)]
    end

    FileUtils.mkdir_p files_dir
    stream = Misc.save_stream(file(:mutation_annotations), mutation_annotations.stream)

    mi_annotations = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => ["Residue"].concat(fields), :type => :double, :namespace => organism
    mi_annotations.init
    TSV.traverse stream, :into => mi_annotations, :bar => "Sorting by Mut. Iso." do |mutation, values|
      mutation = mutation.first if Array === mutation
      Misc.zip_fields(values).each do |mi, *rest|
        mi_annotations.add mi, rest.collect{|v| v.nil? ? [] : v.to_s.split(";") }
      end
      nil
    end
    TSV.collapse_stream mi_annotations.stream
  end

  dep Sequence, :mutated_isoforms_fast
  task :variant_neighbours => :tsv do 
    mutated_isoforms = step(:mutated_isoforms_fast)
    mutated_isoforms.grace
    organism = mutated_isoforms.info[:inputs][:organism]

    neighbours = TSV::Dumper.new :key_field => "Genomic Mutations", :fields => ["Mutated Isoform", "Residue", "PDB", "Neighbours"], :type => :double, :namespace => organism
    neighbours.init
    TSV.traverse mutated_isoforms, :cpus => $cpus, :bar => "Variant neighbours", :into => neighbours do |mutation,mis|
      next if mis.nil? or mis.empty?
      all_mutation_annotations = []
      mis.each do |mi|
        next unless mi =~ /(.*):([A-Z])(\d+)([A-Z])$/
        next if $2 == $4
        isoform = $1
        residue = $3.to_i

        n =  Persist.persist("Neighbours", :yaml, :dir => NEIGHBOURS, :other => {:isoform => isoform, :sequence => residue, :organism => organism}) do 
          Structure.neighbours_i3d(isoform, [residue], organism)
        end
        next if n.empty?
        n.each do |r,v|
          _p,_pos,pdb,n = v
          all_mutation_annotations << [mi,residue,pdb,n]
        end
      end
      annots = Misc.zip_fields(all_mutation_annotations)
      [mutation, annots]
    end

    FileUtils.mkdir_p files_dir
    stream = Misc.save_stream(file(:mutation_neighbours), neighbours.stream)

    mi_annotations = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => ["Residue", "PDB", "Neighbours"], :type => :double
    mi_annotations.init
    TSV.traverse stream, :into => mi_annotations do |mutation, values|
      mutation = mutation.first if Array === mutation
      Misc.zip_fields(values).each do |mi, *rest|
        mi_annotations.add mi, rest.collect{|v| v.nil? ? [] : v.to_s.split(";") }
      end
      nil
    end
    CMD.cmd('sort -u', :in => mi_annotations.stream, :pipe => true)
  end

  #{{{ ANNOTATE NEIGHBOURS
  dep :variant_neighbours
  input :database, :select, "Database of annotations", "UniProt", :select_options => ["UniProt", "COSMIC", "Appris", "InterPro", "variants"]
  task :annotate_neighbours => :tsv do |database|

    variant_neighbours = step(:variant_neighbours)
    mutated_isoforms = variant_neighbours.step(:mutated_isoforms_fast)
    mutated_isoforms.grace
    organism = mutated_isoforms.info[:inputs][:organism]

    case database
    when "COSMIC"
      cosmic_residue_mutations = Structure.COSMIC_residues
      cosmic_mutation_annotations = Structure.COSMIC_mutation_annotations
      fields = ["Genomic Mutation"] + cosmic_mutation_annotations.fields
    when "variants"
      iso2uni = Organism.protein_identifiers(organism).index :target => "UniProt/SwissProt Accession", :persist => true
      iso2sequence = Organism.protein_sequence(organism).tsv :type => :single, :persist => true
      annotations = Structure.UniProt_mutation_annotations
      fields = ["UniProt Variant ID", "SNP ID",  "Type of Variant", "Disease"]
    when "UniProt", "variants"
      iso2uni = Organism.protein_identifiers(organism).index :target => "UniProt/SwissProt Accession", :persist => true
      iso2sequence = Organism.protein_sequence(organism).tsv :type => :single, :persist => true
      fields = ["UniProt Features", "UniProt Feature locations", "UniProt Feature Descriptions"]
    when "InterPro"
      iso2uni = Organism.protein_identifiers(organism).index :target => "UniProt/SwissProt Accession", :persist => true
      iso2sequence = Organism.protein_sequence(organism).tsv :type => :single, :persist => true
      fields = ["InterPro ID", "Range"]
    when "Appris"
      fields = ["Appris Features", "Appris Feature locations", "Appris Feature Descriptions"]
    end

    mutation_annotations = TSV::Dumper.new :key_field => "Genomic Mutations", :fields => ["Mutated Isoform", "Residue"].concat(fields), :type => :double
    mutation_annotations.init
    variant_neighbours.join
    TSV.traverse variant_neighbours.file(:mutation_neighbours), :cpus => $cpus, :bar => "Annot. #{ database }", :into => mutation_annotations do |mutation,values|
      mutation = mutation.first if Array === mutation
      mis, residues, pdb, neighbours = values
      next if mis.nil? or mis.empty? 
      all_mutation_annotations = []
      mis.each do |mi|
        next unless mi =~ /(.*):([A-Z])(\d+)([A-Z])$/
        next if $2 == $4
        isoform = $1
        residue = $3.to_i
        
        n = neighbours ? 
          neighbours.collect{|e| e.split(";") } : 
          [residue-1,residue+1]


        n.flatten.each do |residue|
          residue = residue.to_i

          annotations = case database
                        when "COSMIC"
                          isoform_residue = isoform + ":" << residue.to_s
                          mutations = cosmic_residue_mutations[isoform_residue]
                          next if mutations.nil? 
                          mutations.uniq!
                          next if mutations.empty?
                          annot = cosmic_mutation_annotations.values_at(*mutations)
                          Misc.zip_fields(annot.compact).collect{|v| v * "|" }
                          annot.unshift mutations * "|"

                        when "Appris"
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

                        when "InterPro"
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

                        when "UniProt", "variants"
                          uniprot = iso2uni[isoform]
                          next if uniprot.nil?
                          sequence = iso2sequence[isoform]
                          next if sequence.nil?

                          features =  Persist.persist("Corrected UniProt features", :yaml, :dir => CORRECTED_FEATURES, :other => {:uniprot => uniprot, :sequence => sequence}) do 
                            Structure.corrected_uniprot_features(uniprot, sequence)
                          end
                          next if features.empty?

                          if database == "variants"
                            overlapping = [[],[],[],[]]
                          else
                            overlapping = [[],[],[]]
                          end
                          features.select{|info|
                            if database == "variants"
                              info[:type] == "VARIANT"
                            else
                              case info[:type]
                              when "VAR_SEQ", "CONFLICT", "CHAIN", "UNSURE"
                                false
                              when "DISULFID", "CROSSLNK", "VARIANT"
                                info[:start] == residue or info[:end] == residue
                              else
                                info[:start] <= residue and info[:end] >= residue
                              end
                            end
                          }.each{|info|
                            if database == "variants"
                              if info[:description].match(/(VAR_\d+)/)
                                id = $1
                                next unless annotations.include? id
                                overlapping[0] << id
                                annots = annotations[id]
                                overlapping[1] << annots[2]
                                overlapping[2] << annots[1]
                                overlapping[3] << annots[3]
                              end
                            else
                              overlapping[0] << info[:code]
                              overlapping[1] << [info[:start], info[:end]] * ":"
                              overlapping[2] << (info[:description] || "").strip.sub(/\.$/,'')
                            end
                          }

                          next if overlapping.first.empty?

                          overlapping
                        end

          annotations.unshift residue
          annotations.unshift mi
          all_mutation_annotations << annotations.collect{|v| Array === v ? v * ";" : v }
        end
      end

      next if all_mutation_annotations.empty?
      [mutation, Misc.zip_fields(all_mutation_annotations.uniq)]
    end

    FileUtils.mkdir_p files_dir
    stream = Misc.save_stream(file(:mutation_annotations), mutation_annotations.stream)

    mi_annotations = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => ["Residue"].concat(fields), :type => :double
    mi_annotations.init
    TSV.traverse stream, :into => mi_annotations do |mutation, values|
      mutation = mutation.first if Array === mutation
      Misc.zip_fields(values).each do |mi, *rest|
        mi_annotations.add mi, rest.collect{|v| v.nil? ? [] : v.to_s.split(";") }
      end
      nil
    end
    TSV.collapse_stream mi_annotations.stream
  end


  dep Sequence, :mutated_isoforms_fast
  task :variant_interfaces => :tsv do 
    mutated_isoforms = step(:mutated_isoforms_fast)
    mutated_isoforms.grace
    organism = mutated_isoforms.info[:inputs][:organism]

    interfaces = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => ["Mutated Isoform", "Residue", "Partner Ensembl Protein ID", "PDB", "Partner Residues"], :type => :double
    interfaces.init
    TSV.traverse mutated_isoforms, :cpus => $cpus, :bar => "Variant interfaces", :into => interfaces do |mutation,mis|
      next if mis.nil? or mis.empty?
      all_mutation_annotations = []
      mis.each do |mi|
        next unless mi =~ /(.*):([A-Z])(\d+)([A-Z])$/
        next if $2 == $4
        isoform = $1
        residue = $3.to_i

        n = Persist.persist("Interface neighbours", :yaml, :dir => INTERFACE_NEIGHBOURS, :other => {:isoform => isoform, :sequence => residue, :organism => organism}) do 
          Structure.interface_neighbours_i3d(isoform.dup, [residue], organism)
        end

        next if n.nil? or n.empty?

        n.each do |r,v|
          _pos,part,pdb,n = v
          all_mutation_annotations << [mi,residue,part.uniq,pdb*";",n*";"]
        end
      end
      annots = Misc.zip_fields(all_mutation_annotations.uniq)
      next if annots.empty?
      [mutation, annots]
    end

    FileUtils.mkdir_p files_dir
    stream = Misc.save_stream(file(:mutation_interfaces), interfaces.stream)

    mi_annotations = TSV::Dumper.new :key_field => "Mutated Isoform", :fields => ["Residue", "Partner Ensembl Protein ID", "PDB", "Partner Residues"], :type => :double
    mi_annotations.init
    TSV.traverse stream, :into => mi_annotations do |mutation, values|
      mutation = mutation.first if Array === mutation
      Misc.zip_fields(values).each do |mi, *rest|
        mi_annotations.add mi, rest.collect{|v| v.nil? ? [] : v.to_s.split(";") }
      end
      nil
    end
    CMD.cmd('sort -u', :in => mi_annotations.stream, :pipe => true)
  end
end
