require 'rbbt-util'
module PDBHelper
  def self.pdb_stream(pdb = nil, pdbfile = nil)
    return StringIO.new(pdbfile) if (pdb.nil? or pdb.empty?) and not pdbfile.nil? and not pdbfile.empty?
    return Open.open(pdb) if pdb and (Open.remote?(pdb) or Open.exists?(pdb))
    return Open.open("http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=#{pdb}") unless pdb.nil?

    raise "No valid pdb provided: #{ pdb }"
  end

  def self.atoms(pdb = nil, pdbfile = nil)
    io = pdb_stream(pdb,pdbfile)
    str = ""
    while line = io.gets and not line =~ /^END/
      str << line if line =~ /^ATOM/
    end
    str
  end

  def self.pdb_chain_sequences(pdb = nil, pdbfile = nil)
    atoms = PDBHelper.atoms(pdb, pdbfile)

    chains = {}
    atoms.split("\n").each do |line|
      chain = line[20..21].strip
      aapos = line[22..25].to_i
      aa    = line[17..19]

      next if aapos < 0

      chains[chain] ||= Array.new
      chains[chain][aapos-1] = aa
    end

    chains.each do |chain,chars|
      chains[chain] = chars.collect{|aa| aa.nil? ? '?' : Misc::THREE_TO_ONE_AA_CODE[aa.downcase]} * ""
    end

    chains
  end

  def self.pdb_atom_distance(distance, pdb = nil, pdbfile = nil)
    atoms = self.atoms(pdb, pdbfile)
    atom_positions = {}
    atoms.split("\n").each{|line|
      code = line[13..26]
      x = line[30..37].to_f
      y = line[38..45].to_f
      z = line[46..53].to_f

      atom_positions[code] = [x,y,z]
    }
    atom_positions
    atoms = atom_positions.keys.sort_by{|atom| atom[9..13].to_i}

    atom_distances = []
    atoms.each_with_index do |atom1,i|
      position1 = atom_positions[atom1]
      atoms[i+1..-1].each do |atom2|
        next if (atom1[9..13] == atom2[9..13]) 
        next if ((atom1[9..13].to_i == atom2[9..13].to_i + 1) and 
        ((atom1[0] == "C" and atom2[0] == "N") or 
          (atom1[0] == "0" and atom2[0] == "N") or 
          (atom1[0] == "C" and atom2[0..1] == "CA")))

        position2 = atom_positions[atom2]
        dx = position1[0] - position2[0]
        dy = position1[1] - position2[1]
        dz = position1[2] - position2[2]
        next if dx.abs > 5 or dy.abs > 5 or dz.abs > 5
        dist = Math.sqrt(dx**2 + dy**2 + dz**2)
        next if dist > 5
        atom_distances << [atom1, atom2, dist]
      end
    end
    atom_distances
  end

  def self.pdb_close_residues(distance, pdb = nil, pdbfile = nil)
    atom_distances = pdb_atom_distance(distance, pdb, pdbfile)
    close_residues = {}
    atom_distances.each do |atom1, atom2, dist|
      aa1 = atom1.split(/\s+/).values_at(2,3) * ":"
      aa2 = atom2.split(/\s+/).values_at(2,3) * ":"
      close_residues[aa1] ||= []
      close_residues[aa1] << aa2
    end

    close_residues.each do |aa1, list| list.uniq! end
    close_residues
  end
end
