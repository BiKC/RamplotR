const assert = require("assert");
const mapping = require("../shinyRam/www/canonical-mapping.js");

const payload = {
  "1abc": {
    UniProt: {
      P12345: {
        identifier: "TEST_HUMAN",
        mappings: [{
          entity_id: 1,
          chain_id: "A",
          struct_asym_id: "A",
          unp_start: 10,
          unp_end: 12,
          identity: 1,
          coverage: 0.5,
          start: {
            residue_number: 1,
            author_residue_number: 101,
            author_insertion_code: ""
          },
          end: {
            residue_number: 3,
            author_residue_number: 103,
            author_insertion_code: ""
          }
        }]
      }
    }
  }
};

const rows = mapping.normalizePayload(payload, "1ABC");
assert.strictEqual(rows.length, 1);
assert.deepStrictEqual(rows[0], {
  pdb_id: "1ABC",
  uniprot_accession: "P12345",
  uniprot_identifier: "TEST_HUMAN",
  entity_id: 1,
  chain: "A",
  struct_asym_id: "A",
  unp_start: 10,
  unp_end: 12,
  pdb_start: 1,
  pdb_end: 3,
  author_start: 101,
  author_end: 103,
  author_start_insertion: "",
  author_end_insertion: "",
  identity: 1,
  coverage: 0.5
});

assert.strictEqual(mapping.normalizePayload({}, "1ABC").length, 0);
assert.strictEqual(mapping.entryFromPayload({"1ABC": {x: 1}}, "1abc").x, 1);

console.log("Canonical mapping browser normalization tests passed.");
