module Structure
  def self.uni2iso(organism = Organism.default_code("Hsa"))
    @@uni2iso ||= {}
    #@@uni2iso[organism] ||= Organism.protein_identifiers(organism).index :fields => ["UniProt/SwissProt Accession"], :target => "Ensembl Protein ID", :persist => true, :unnamed => true
    @@uni2iso[organism] ||= Organism.uniprot2ensembl(organism).index :fields => ["UniProt/SwissProt Accession"], :target => "Ensembl Protein ID", :persist => true, :unnamed => true
  end

  def self.iso2uni(organism = Organism.default_code("Hsa"))
    @@iso2uni ||= {}
    #@@iso2uni[organism] ||= Organism.protein_identifiers(organism).index :target => "UniProt/SwissProt Accession", :fields => ["Ensembl Protein ID"], :persist => true, :unnamed => true
    @@iso2uni[organism] ||= Organism.ensembl2uniprot(organism).index :target => "UniProt/SwissProt Accession", :fields => ["Ensembl Protein ID"], :persist => true, :unnamed => true
  end

  def self.iso2seq(organism = Organism.default_code("Hsa"))
    @@iso2seq ||= {}
    @@iso2seq[organism] ||= Organism.protein_sequence(organism).tsv :persist => true, :unnamed => true
  end

end
