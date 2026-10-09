const assert=require("node:assert/strict");
require("../shinyRam/www/atlas-sifts-exact.js");
const ligand=require("../shinyRam/www/atlas-ligand-contacts.js");
const exact={rows:[
  {observed:true,struct_asym_id:"X",label_seq_id:1,uniprot_resi:101},
  {observed:true,struct_asym_id:"X",label_seq_id:2,uniprot_resi:102}
]};
const atomHeader=[
"group_PDB","label_atom_id","label_asym_id","label_seq_id",
"label_comp_id","type_symbol","Cartn_x","Cartn_y","Cartn_z",
"pdbx_PDB_model_num","label_alt_id","auth_seq_id","pdbx_PDB_ins_code"];
const rows=[
"ATOM N X 1 ALA N 0 0 0 1 . 101 .",
"ATOM CA X 1 ALA C 1 0 0 1 . 101 .",
"ATOM CA X 2 GLY C 30 0 0 1 . 102 .",
"HETATM C1 L . ATP C 3 0 0 1 . 900 .",
"HETATM H1 L . ATP H 0 0 0 1 . 900 .",
"HETATM MG I . MG MG 50 0 0 1 . 10 .",
"HETATM O W . HOH O 0 0 0 1 . 1 .",
"HETATM C1 X 2 MSE C 1 0 0 1 . 102 .",
"HETATM C1 Q . GTP C 2 0 0 2 . 9 .",
"HETATM C1 B . FAD C 2 0 0 1 B 10 ."];
const cif="data_test\n#\nloop_\n"+atomHeader.map(h=>"_atom_site."+h).join("\n")+
  "\n"+rows.join("\n")+"\n#";
const out=ligand.extract(cif,exact);
assert.equal(out.status,"measured");
assert.equal(out.total_nonwater_sites,2);
assert.equal(out.nearby_sites,1);
assert.equal(out.sites[0].comp_id,"ATP");
assert.equal(out.sites[0].nearest_uniprot_resi,101);
assert.equal(out.sites[0].min_distance_A,2);
assert.equal(out.sites[0].heavy_atoms,1);
assert.equal(out.radius_A,4.5);
const noLigand=ligand.extract(cif.replace(/HETATM (?:C1 L|H1 L|MG I).*\n/g,""),exact);
assert.equal(noLigand.status,"measured");
assert.equal(noLigand.nearby_sites,0);
assert.equal(noLigand.total_nonwater_sites,0);
assert.match(noLigand.warning,/does not establish an apo state/);
assert.equal(ligand.extract(cif,{rows:[]}).status,"unavailable");
assert.throws(()=>ligand.extract(cif.replace("_atom_site.type_symbol",
  "_atom_site.unknown"),exact),/lacks fields/);
assert.throws(()=>ligand.extract(cif,exact,20),/Invalid/);
console.log("Exact-mapped deposited hetero heavy-atom contact tests passed.");
