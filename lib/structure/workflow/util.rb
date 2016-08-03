module Structure

  def self.name2pi (mutation, organism)
   @@index = Organism.identifiers(organism).index :target => "Ensembl Gene ID", :order => true, :persist => true
   @@ensg2enst = Organism.transcripts(organism).tsv :key_field => "Ensembl Gene ID", :fields => ["Ensembl Transcript ID"], :persist => true, :merge => true, :type => :flat
   @@enst2ensp = Organism.transcripts(organism).index :target => "Ensembl Protein ID", :fields => ["Ensembl Transcript ID"], :persist => true, :merge => true
   @@enst2name = Organism.transcript_name(organism).index :target => "Ensembl Transcript Name", :fields => ["Ensembl Transcript ID"], :persist => true, :merge => true
   @@uni2ensp = Organism.uniprot2ensembl(organism).index :persist => true, :target => "Ensembl Protein ID"
   @@ensp2uni = Organism.ensembl2uniprot(organism).index :persist => true, :target => "UniProt/SwissProt Accession"

  orig_gene, _sep, change = mutation.partition ":"
  gene = @@index[orig_gene]

  return nil if gene.nil? or @@ensg2enst[gene].nil?

  gene_transcripts = @@ensg2enst[gene].sort_by{|t| @@enst2name[t].split("-").last.to_i}
  gene_isoforms = @@enst2ensp.values_at(*gene_transcripts).compact.reject{|i| i.empty?}
  perfect_isoforms = gene_isoforms.select{|p|
    _uni = @@ensp2uni[p]
    _p = @@uni2ensp[_uni]
    p == _p
  }

  principal_transcripts = (Appris::PRINCIPAL_TRANSCRIPTS & gene_transcripts)
  principal_isoforms = @@enst2ensp.values_at *principal_transcripts
  uni_pricipal_isoforms = principal_isoforms.select{|p| @@ensp2uni[p]}

  perfect_principal_isoforms = uni_pricipal_isoforms & perfect_isoforms

  if perfect_principal_isoforms.any?
    protein = perfect_principal_isoforms.first
  elsif uni_pricipal_isoforms.any?
    protein = uni_pricipal_isoforms
  elsif perfect_isoforms.any?
    protein = perfect_isoforms.first
  elsif principal_isoforms.any?
    protein = principal_isoforms.first
  else
    protein = gene_isoforms.first
  end

  return nil if protein.nil? or protein.empty?

  [protein, change] * ":"
  end
end

