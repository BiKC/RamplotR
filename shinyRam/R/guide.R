# In-app, task-oriented help. Keep this in sync with UI controls and the
# scientific limitations documented in docs/. No new runtime packages.
ram_guide_ui <- function() {
  step <- function(number, title, description) {
    tags$li(
      tags$span(class = "ram-guide-step-number", as.character(number)),
      tags$div(tags$strong(title), tags$p(description))
    )
  }
  jump <- function(id, label) {
    actionLink(id, label, class = "ram-guide-go")
  }
  section <- function(id, eyebrow, title, intro, steps, note = NULL, go = NULL) {
    tags$section(id = id, class = "ram-guide-section",
      tags$div(class = "ram-guide-section-head",
        tags$span(class = "ram-guide-eyebrow", eyebrow),
        tags$h3(title),
        tags$p(intro)
      ),
      tags$ol(class = "ram-guide-steps", steps),
      if (!is.null(note))
        tags$div(class = "ram-guide-note", note),
      if (!is.null(go)) go
    )
  }

  tags$div(class = "ram-subtab-content ram-guide",
    tags$div(class = "ram-guide-intro",
      tags$p(class = "ram-guide-eyebrow", "Documentation and tutorials"),
      tags$h2("What do you want to find out?"),
      tags$p(paste(
        "You do not need to use every feature. Start with a research question,",
        "then follow the relevant steps. Each workflow below refers to controls",
        "that are already in RamplotR."
      )),
      tags$div(class = "ram-guide-quick-start",
        tags$strong("New to RamplotR?"),
        tags$span(paste(
          "Load PDB 1CRN, choose Analyze structure, and click a residue",
          "in the Ramachandran plot. Its details and 3D position will follow",
          "your selection."
        )),
        jump("guideGoPlot", "Open the Ramachandran plot")
      )
    ),
    tags$nav(class = "ram-guide-index", "aria-label" = "Guide sections",
      tags$a(href = "#ram-guide-inspect", "Inspect a residue"),
      tags$a(href = "#ram-guide-prediction", "Predicted structures"),
      tags$a(href = "#ram-guide-pair", "Compare two structures"),
      tags$a(href = "#ram-guide-groups", "Compare groups"),
      tags$a(href = "#ram-guide-atlas", "Experimental Atlas"),
      tags$a(href = "#ram-guide-export", "Export and report"),
      tags$a(href = "#ram-guide-interpretation", "Interpret the evidence")
    ),
    tags$div(class = "ram-guide-overview",
      tags$article(
        tags$span("ONE MODEL"), tags$strong("Is this residue unusual?"),
        tags$p("Inspect native RamplotR density regions, Rama8000 validation and local 3D context.")
      ),
      tags$article(
        tags$span("TWO MODELS"), tags$strong("Where did the conformation change?"),
        tags$p("Align chains, measure wrapped backbone angle differences and inspect residues in both structures.")
      ),
      tags$article(
        tags$span("MODEL SETS"), tags$strong("Do predictions or groups agree?"),
        tags$p("Check per-residue conformational agreement across models, or contrast two named structure sets.")
      ),
      tags$article(
        tags$span("EXPERIMENTAL EVIDENCE"), tags$strong("Are other states observed?"),
        tags$p("Discover PDB entries using UniProt, verify SIFTS mappings and explore experimental geometry groups.")
      )
    ),
    section("ram-guide-inspect", "01 · First analysis",
      "Inspect one protein and its residues",
      "Start here for experimental structures, deposited PDB files, and single predicted models.",
      tagList(
        step(1, "Choose a source",
          "At the top of the page, choose PDB ID, Upload file or AlphaFold DB. Try 1CRN as a small PDB example."),
        step(2, "Click Analyze structure",
          "The Ramachandran plot and 3D viewer use the selected protein chains. Choose the reference dataset and region classification in the settings sidebar."),
        step(3, "Select a residue",
          "Click a Ramachandran point, a sequence-navigator position or a residue-list row. The selection is linked to the 3D structure and residue inspector."),
        step(4, "Inspect the evidence",
          "Check φ/ψ, native RamplotR region, independent six-class Rama8000 category, and available local geometry or nearby ligands. The issue queue helps focus on unusual residues.")
      ),
      tags$p("The native RamplotR density label and Rama8000 Favored/Allowed/Outlier are different classifications. They are intentionally shown separately."),
      jump("guideGoResidues", "Open the residue list")
    ),
    section("ram-guide-prediction", "02 · Model confidence and ensembles",
      "Investigate AlphaFold, ColabFold, AF3 or ESMFold",
      "Use this when your question is whether multiple predictions agree on a local backbone conformation.",
      tagList(
        step(1, "Declare the correct prediction source",
          "For an uploaded model, open Prediction settings near the upload form and choose AF2/ColabFold, AF3, ESMFold or another predicted model. AlphaFold DB imports are identified automatically."),
        step(2, "Add confidence data when available",
          "Load optional compatible PAE or AF3 confidence JSON files in Prediction settings. Per-residue pLDDT and pairwise PAE are different measures."),
        step(3, "Inspect a single prediction",
          "After analysis, use the linked residue inspector and confidence views. A high pLDDT does not establish that the protein adopts a single biological state."),
        step(4, "Explore several models",
          "In the Summary tab, expand Ensemble analysis. Use a multi-model structure or the Prediction ensemble import to compare models, seeds or samples where supported."),
        step(5, "Check disagreement",
          "Look for residues with circular φ/ψ dispersion, alternative backbone states or differences between pLDDT and model-to-model agreement.")
      ),
      tags$p("Prediction samples are not independent thermodynamic observations. Missing PAE is not zero uncertainty, and pLDDT does not validate a conformational transition."),
      jump("guideGoSummary", "Open Summary and ensemble analysis")
    ),
    section("ram-guide-pair", "03 · Conformational Change Explorer",
      "Compare two structures residue by residue",
      "Useful for apo/holo, mutant/wild type, experimental/predicted, and different modeling methods.",
      tagList(
        step(1, "Load the first structure",
          "It becomes your initial reference in the main input form."),
        step(2, "Open Compare",
          "Expand Comparison structure, choose a second PDB accession or upload a PDB/mmCIF file, and load it."),
        step(3, "Review chain alignment",
          "RamplotR aligns chains by sequence and reports identity and reference coverage. Residue numbers are not assumed to match."),
        step(4, "Inspect local changes",
          "Use the aligned residue track, the angle plot and the comparison table. Click a residue to highlight both corresponding positions in 3D."),
        step(5, "Swap or filter",
          "Swap primary/comparison roles when useful. Filters include wrapped Δφ/Δψ, broad backbone-state transitions, validation changes and prediction confidence where available.")
      ),
      tags$p("Δφ and Δψ wrap around ±180°. A large angular change is a structural observation, not proof of a functional switch."),
      jump("guideGoCompare", "Open Compare")
    ),
    section("ram-guide-groups", "04 · Group conformational comparison",
      "Compare two sets of related protein structures",
      "Use this when several structures belong to each condition, not when you only have one pair.",
      tagList(
        step(1, "Load a representative reference",
          "Load a protein first, then open Compare and expand Compare groups of structures. Pick its Reference chain."),
        step(2, "Set alignment requirements",
          "The default minimum chain identity and reference coverage are both 70%. These exclude poorly matched uploaded chains; they are not biological-state thresholds."),
        step(3, "Build Group A and Group B",
          "The loaded structure can count as a Group A member. Upload additional Group A structures if needed, and upload at least one Group B structure. Label each condition clearly."),
        step(4, "Click Analyse groups",
          "Each uploaded file contributes model 1, with its best-matching protein chain. RamplotR reports which uploaded structures were matched or omitted."),
        step(5, "Review support, dispersion and exports",
          "Inspect group-level circular φ/ψ summaries, changed-residue evidence and the matched-members table. Export residue-level and matched-chain CSVs after a successful run.")
      ),
      tags$p("Alternatively, use verified Atlas groups directly: after calculating experimental geometry, select disjoint Group A and Group B members in Atlas and send them to Compare Groups. The exact observed UniProt backbone cache is reused without uploads. Classification metrics are unavailable in this mode. A single member per condition is a comparison, not a distribution; treat support and dispersion cautiously."),
      jump("guideGoGroups", "Open group comparison")
    ),
    section("ram-guide-atlas", "05 · Experimental Conformational Atlas",
      "Discover and compare experimental structures",
      "Use the Atlas to find experimental counterparts for a UniProt protein and explore possible conformational differences.",
      tagList(
        step(1, "Verify archive connectivity if needed",
          "Use Check browser access to RCSB and PDBe, then Test archive connections. This runs directly from the current browser and checks RCSB Search, polymer metadata, and PDBe updated mmCIF/SIFTS on a known protein. Failures are shown independently; the diagnostic does not change your analysis."),
        step(2, "Search a UniProt accession",
          "In Atlas, select Find experimental structures. RCSB returns experimental PDB polymer entities in pages, with metadata such as method and resolution."),
        step(3, "Verify exact residue mapping",
          "Click Verify SIFTS mapping on candidate entries. RamplotR retrieves PDBe updated mmCIF and checks exact PDB-to-UniProt residue correspondence, including insertion codes."),
        step(4, "Build a verified cohort",
          "Once entries are verified, review distinct observed UniProt positions, incomplete mappings, conflicts and downloadable residue coverage."),
        step(5, "Review experimental construct compatibility",
          "For the selected structures, inspect pairwise observed UniProt core coverage and the monomer identities encoded in the mmCIF polymer scheme. Missing or unmatched positions do not establish equivalence. Any confirmed residue-chemistry differences require your explicit acknowledgement before exploratory grouping; export the review CSV for provenance."),
        step(6, "Let Atlas suggest geometric clusters",
          "Select at least two exact-verified experimental entities, review construct chemistry and choose Cluster experimental structures. Automatic mode searches for coherent groupings using Cα distance-map separation and mean silhouette, and can recommend no split when evidence is weak. Its numeric checks are exploratory, not calibrated functional-state labels."),
        step(7, "Review or adjust cluster membership",
          "Examine the dendrogram, group sizes and the Why did Atlas suggest these groups? quality details. You can switch to a manual distance cutoff or edit the suggested Group A/Group B memberships. Use these entries in Compare Groups passes verified experimental backbone torsions directly without file uploads."),
        step(8, "Investigate local differences",
          "Choose two geometric-group representatives, then select Find local backbone changes. Click a canonical residue to see paired φ/ψ and exact PDB residue identifiers."),
        step(9, "Inspect both representative structures in 2D/3D",
          "Click Inspect both structures in 2D/3D. Compare now uses shared, observed UniProt positions from the exact verified SIFTS mappings, even when author numbering or constructs differ. Unmapped residues are excluded and the compared coverage is shown. A missing or ambiguous mapping is never replaced by a guessed match.")
      ),
      tags$p("Atlas groups are exploratory geometric similarities, not confirmed functional states. The new chemistry review detects mapped differences and unknowns, but cannot prove isoform/construct equivalence or rule out ligand, experimental-condition or crystal-packing effects. Automatic silhouette and distance-gap safeguards, the optional manual 1.5 Å cutoff and the 30° local-change criterion are not biologically validated state boundaries."),
      jump("guideGoAtlas", "Open Atlas")
    ),
    section("ram-guide-export", "06 · Figures and reproducible output",
      "Save what you can defend",
      "Keep structure provenance, chain mapping and completeness alongside your figures.",
      tagList(
        step(1, "Export a figure",
          "In Summary, expand Export figures and a reproducible report to download a vector SVG or high-resolution PNG using your selected palette."),
        step(2, "Export residue evidence",
          "Residue list, Compare, group comparisons and Atlas provide CSV exports appropriate to their analysis. Exports describe the current matched or filtered data."),
        step(3, "Save a standalone report",
          "The Summary HTML report captures the residue analysis and available provenance."),
        step(4, "Run larger batches offline",
          "For directory-level processing, use scripts/ramplotr-batch.R from the repository; the browser UI is designed for interactive inspection.")
      ),
      tags$p("Always report experimental/prediction provenance, chain identity, usable residue coverage and the distinction between native RamplotR and Rama8000 classifications."),
      jump("guideGoSummary2", "Open export options")
    ),
    tags$section(id="ram-guide-interpretation",class="ram-guide-section",
      tags$span(class="ram-guide-eyebrow", "Reference"),
      tags$h3("Reading the evidence correctly"),
      tags$div(class="ram-guide-faq",
        tags$details(
          tags$summary("Why do the native RamplotR region and Rama8000 category differ?"),
          tags$p("They use different residue classes, reference distributions and categories. RamplotR's native density regions are exploratory; six-class Rama8000 follows the standardized validation model. Do not translate labels one-to-one.")
        ),
        tags$details(
          tags$summary("Does a 30° change mean a real functional switch?"),
          tags$p("No. It flags a candidate for inspection. Same-state experimental controls can also contain changes above 30°, and some domain motions involve relatively little local torsion change.")
        ),
        tags$details(
          tags$summary("Why are some residues missing from comparison results?"),
          tags$p("Structures can have gaps, missing atoms, incomplete peptide geometry, different constructs or uncertain alignment. The app does not invent φ/ψ angles or UniProt positions.")
        ),
        tags$details(
          tags$summary("Are experimental geometric groups validated conformational states?"),
          tags$p("Not yet. They depend on canonical overlap and a user-selected distance threshold. Biological context, construct equivalence and independent experimental controls still need inspection.")
        ),
        tags$details(
          tags$summary("How should I read pLDDT, PAE and prediction-ensemble spread?"),
          tags$p("pLDDT measures local predicted confidence, PAE estimates relative-position uncertainty, and ensemble angle spread measures agreement between sampled predictions. None is interchangeable with experimental validation or a thermodynamic state population.")
        )
      ),
      tags$p(class="ram-guide-more",
        "For algorithms, input formats, deployment and benchmark protocols, see the ",
        tags$a(href="https://github.com/BiKC/RamplotR",
          target="_blank",rel="noopener noreferrer","project documentation"),
        "."
      )
    )
  )
}
