require 'rbbt-util'
module PDBHelper

  PDB_FILES = Structure.cache_dir.pdb_files.find
  Open.repository_dirs << PDB_FILES unless Open.repository_dirs.include? PDB_FILES
  def self.pdb_stream(pdb = nil, pdbfile = nil)
    return StringIO.new(pdbfile) if (pdb.nil? or pdb.empty?) and not pdbfile.nil? and not pdbfile.empty?
    raise "No PDB specified" if pdb.nil? or pdb.empty?
    if Open.exists?(pdb)
      Open.open(pdb)
    else
      begin
        content = Persist.persist("PDB file", :marshal, :persist => true, :dir => PDB_FILES, :other => {:pdb => pdb}) do 
          if pdb.length > 5 and Open.remote?(pdb)
            Open.read(pdb)
          else
            pdb = pdb.sub(/^=/,'')
            Open.read("http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=#{pdb}")
          end
        end
        StringIO.new(content)
      rescue Exception
        Log.exception $!
        raise "No valid pdb provided: #{ Misc.fingerprint pdb }"
      end
    end
  end

  #def self.pdb_stream(pdb = nil, pdbfile = nil)
  #  begin
  #    return StringIO.new(pdbfile) if (pdb.nil? or pdb.empty?) and not pdbfile.nil? and not pdbfile.empty?
  #    return Open.open(pdb) if pdb and (Open.remote?(pdb) or Open.exists?(pdb))
  #    pdb = pdb.sub(/^=/,'')
  #    return Open.open("http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=#{pdb}") unless pdb.nil?

  #  rescue Exception
  #    raise "No valid pdb provided: #{ Misc.fingerprint pdb }"
  #  end

  #  raise "No valid pdb provided: #{ Misc.fingerprint pdb }"
  #end

  CHAIN_SEQUENCES = Structure.cache_dir.chain_sequences.find
  Open.repository_dirs << CHAIN_SEQUENCES unless Open.repository_dirs.include? CHAIN_SEQUENCES
  def self.pdb_chain_sequences(pdb = nil, pdbfile = nil)
    Misc.insist do
      Persist.persist("Chain sequences", :marshal, :persist => true, :dir => CHAIN_SEQUENCES, :other => {:pdb => pdb, :pdbfile => pdbfile}) do 
        chains = {}
        stream = pdb_stream(pdb, pdbfile)
        while line = stream.gets
          break if line =~ /^END/
            next unless line =~ /^ATOM/ and line.length > 21
            chain = line[20..21].strip
            aapos = line[22..25].to_i
            aa    = line[17..19]

            next if aapos <= 0

            chain = "A" if chain == ""
            chains[chain] ||= Array.new
            chains[chain][aapos-1] = aa
        end
        stream.close

        chains.each do |chain,chars|
          chains[chain] = chars.collect{|aa| aa.nil? ? '?' : Misc::THREE_TO_ONE_AA_CODE[aa.downcase]} * ""
        end

        chains
      end
    end
  end

  def self.pdb_atom_distance(distance, pdb = nil, pdbfile = nil)
    stream = pdb_stream(pdb, pdbfile)

    atom_positions = {}
    while line = stream.gets
      break if line =~ /^END/
        #next unless line =~ /^ATOM/
        next unless line =~ /^\w*ATO?M/ and line.length > 53
        code = line[13..26]
        x = line[30..37].to_f
        y = line[38..45].to_f
        z = line[46..53].to_f
        num = code[9..13].to_i

      atom_positions[code] = [x,y,z,num]
    end

    atom_positions

    atoms = atom_positions.
      sort_by{|a,values| values[3] }.
      collect{|atom,v| atom }

    atom_distances = []
    while atom1 = atoms.shift
      position1 = atom_positions[atom1]

      atoms.each do |atom2|
        position2 = atom_positions[atom2]
        next if (position1[3] == position2[3]) 
        next if ((position1[3] == position2[3] + 1) and 
        ((atom1[0] == "C" and atom2[0] == "N") or 
          (atom1[0] == "0" and atom2[0] == "N") or 
          (atom1[0] == "C" and atom2[0..1] == "CA")))

        position2 = atom_positions[atom2]
        dx = position1[0] - position2[0]
        dy = position1[1] - position2[1]
        dz = position1[2] - position2[2]

        next if dx.abs > distance or dy.abs > distance or dz.abs > distance
        dist = Math.sqrt(dx**2 + dy**2 + dz**2)
        next if dist > distance
        atom_distances << [atom1, atom2, dist]
      end

    end

    atom_distances
  end

  def self.pdb_close_residues(distance, pdb = nil, pdbfile = nil)
    Log.low "Computing atom distances (#{ distance }): #{pdb || "pdbfile"}"
    atom_distances = pdb_atom_distance(distance, pdb, pdbfile)

    close_residues = {}
    Log.low "Computing residue distances (#{distance}): #{pdb || "pdbfile"}"
    atom_distances.each do |atom1, atom2, dist|
      aa1 = atom1.strip.split(/\s+/).values_at(2,3) * ":"
      aa2 = atom2.strip.split(/\s+/).values_at(2,3) * ":"
      aa1 = "A" + aa1 if aa1 =~ /^\d/
      aa2 = "A" + aa2 if aa2 =~ /^\d/
      aa1 = aa1[0] + ":" + aa1[1..-2] if aa1[-1] == ":"
      aa2 = aa2[0] + ":" + aa2[1..-2] if aa2[-1] == ":"
      close_residues[aa1] ||= []
      close_residues[aa1] << aa2
    end

    close_residues.each do |aa1, list| list.uniq! end

    close_residues
  end
end
