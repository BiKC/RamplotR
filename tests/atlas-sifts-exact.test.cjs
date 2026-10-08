const assert = require("node:assert/strict");
const sifts = require("../shinyRam/www/atlas-sifts-exact.js");
const mmcif = `data_1abc
#
loop_
_struct.title
'Example with quoted value'
#
loop_
_pdbx_poly_seq_scheme.asym_id
_pdbx_poly_seq_scheme.entity_id
_pdbx_poly_seq_scheme.seq_id
_pdbx_poly_seq_scheme.pdb_strand_id
_pdbx_poly_seq_scheme.pdb_seq_num
_pdbx_poly_seq_scheme.auth_seq_num
_pdbx_poly_seq_scheme.pdb_ins_code
X 1 1 A 101 101 .
X 1 2 A 101 101 A
X 1 3 A 103 103 .
Y 2 1 B 201 201 .
#
loop_
_pdbx_sifts_xref_db.entity_id
_pdbx_sifts_xref_db.asym_id
_pdbx_sifts_xref_db.seq_id
_pdbx_sifts_xref_db.unp_acc
_pdbx_sifts_xref_db.unp_num
_pdbx_sifts_xref_db.observed
1 X 1 P12345 20 1
1 X 2 P12345 21 1
1 X 3 P12345 25 0
2 Y 1 Q99999 8 1
#`;
const result = sifts.exactRows(mmcif,"1abc","1","P12345");
assert.equal(result.matched_sifts_rows,3);
assert.equal(result.unlinked_sifts_rows,0);
assert.equal(result.rows.length,3);
assert.deepEqual(result.rows.map(x=>x.uniprot_resi),[20,21,25]);
assert.deepEqual(result.rows.map(x=>x.resi),[101,101,103]);
assert.deepEqual(result.rows.map(x=>x.insertion_code),["","A",""]);
assert.equal(result.rows[2].observed,false);
assert.equal(result.unobserved_rows,1);
assert.deepEqual(result.rows.map(x=>x.chain),["A","A","A"]);
assert.equal(sifts.exactRows(mmcif,"1abc","2","Q99999").rows[0].chain,"B");
assert.throws(()=>sifts.exactRows(mmcif,"../a","1","P12345"),/Invalid/);
assert.throws(()=>sifts.exactRows(mmcif.replace(
"_pdbx_poly_seq_scheme.pdb_ins_code",
"_pdbx_poly_seq_scheme.unknown"),"1abc","1","P12345"),/lacks/);
const duplicate = mmcif.replace("1 X 2 P12345 21 1",
  "1 X 2 P12345 21 1\n1 X 2 P12345 21 1");
assert.equal(sifts.exactRows(duplicate,"1abc","1","P12345").rows.length,3);
// No interpolated 102: the PDB-author 101A and 103 identifiers are literal.
assert.ok(!result.rows.some(x=>x.resi===102));
assert.deepEqual(sifts.tokenize('  \'A title\'  "other title" # comment'),
  ["A title","other title"]);
// Coordinates are joined by label_asym_id + label_seq_id, never author
// residue numbers. Extra models/altlocs and non-protein atoms are excluded.
const withAtoms = mmcif + `
loop_
_atom_site.group_PDB
_atom_site.label_atom_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
_atom_site.label_alt_id
ATOM CA X 1 1.0 2.0 3.0 1 .
ATOM CA X 2 4.0 5.0 6.0 1 A
ATOM CA X 3 7.0 8.0 9.0 1 .
ATOM CA X 1 99.0 99.0 99.0 2 .
ATOM N X 2 10.0 10.0 10.0 1 .
HETATM CA Y 1 20.0 20.0 20.0 1 .
#`;
const ca = sifts.extractMappedCA(withAtoms,result);
assert.equal(ca.length,2,"Unobserved SIFTS residues must not be sent as C-alpha points.");
assert.deepEqual(ca.map(x=>x.label_seq_id),[1,2]);
assert.deepEqual([ca[0].x,ca[0].y,ca[0].z],[1,2,3]);
assert.equal(result.rows[1].label_seq_id,2);
assert.throws(()=>sifts.extractMappedCA(withAtoms.replace(
  "_atom_site.Cartn_z","_atom_site.unknown"),result),/lacks/);
console.log("Exact SIFTS mmCIF parser tests passed.");
