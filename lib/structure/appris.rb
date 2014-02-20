module Structure
  def self.isoform2transcript
    @isoform2transcript ||= Organism.transcripts("Hsa").index :target => "Ensembl Transcript ID", :fields => "Ensembl Protein ID", :persist => true
  end

  def self.appris_dataset
    @appris_dataset ||= begin
                          Rbbt.share.data["appris_data.fir_spa_thu_cra.gen19.v2.tsv"].find(:lib).tsv :key_field => "Ensembl Transcript ID", :persist => true
                        end
  end

  def self.appris_features(isoform)
    transcript = isoform2transcript[isoform]

    values = appris_dataset[transcript]
    return [] if values.nil?

    features = []

    values["firestar functional residues"].each do |res|
      position, ligand = res.split ":"
      features << {:type => "firestar", :start => position.to_i, :end => position.to_i, :description => ligand}
    end

    values["spade whole domain"].each do |res|
      position, pfam_acc = res.split ":"
      start, eend = position.split("-")
      pfam_acc.sub!(/\.\d+$/, '')
      features << {:type => "spade", :start => start.to_i, :end => eend.to_i, :description => pfam_acc}
    end

    values["thump transmembrane helix"].each do |res|
      position, damage = res.split ":"
      start, eend = position.split("-")
      damage = damage == "1" ? "Damaged" : "Normal"
      features << {:type => "thump", :start => start.to_i, :end => eend.to_i, :description => damage}
    end

    values["crash signal peptide"].each do |res|
      start, eend = res.split("-")
      features << {:type => "crash", :start => start.to_i, :end => eend.to_i, :description => "Signal peptide"}
    end

    features
  end
end
