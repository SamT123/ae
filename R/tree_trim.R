#! /usr/bin/env Rscript

main = function() {
  message("TRIM TREE")

  .libPaths(
    c(
      "/home/samt123/R/x86_64-pc-linux-gnu-library/4.2",
      .libPaths()
    )
  )

  library(seqUtils, quietly = T, warn.conflicts = F, verbose = F)
  library(convergence, quietly = T, warn.conflicts = F, verbose = F)

  library(treeio, quietly = T, warn.conflicts = F, verbose = F, )

  library(purrr, quietly = T, warn.conflicts = F, verbose = F)
  library(dplyr, quietly = T, warn.conflicts = F, verbose = F)
  library(stringr, quietly = T, warn.conflicts = F, verbose = F)

  convergence:::add_to_PATH("/home/sat65/miniforge3/envs/treebuild/bin") # usher

  # args ------------------------------------------------------------
  args <- commandArgs(trailingOnly = TRUE)

  subtype = args[[1]]
  inphy.txt = args[[2]]
  infasta = args[[3]]

  inphy = readLines(inphy.txt, n = 1)

  outphy = paste0(
    paste0(
      fs::path_ext_remove(inphy),
      "-trim."
    ),
    fs::path_ext(inphy)
  )

  outphy.txt = paste0(
    paste0(
      fs::path_ext_remove(inphy.txt),
      "-trim."
    ),
    fs::path_ext(inphy.txt)
  )

  outfasta = paste0(
    paste0(
      fs::path_ext_remove(infasta),
      "-trim."
    ),
    fs::path_ext(infasta)
  )

  if (all(fs::file_exists(outphy, outphy.txt, outfasta))) {
    message("All output files already exist. Not remaking.")
  }

  # specify trim ----------------------------------------------------
  if (subtype == "h3") {
    trim = T
    strain = "COTE_DIVOIRE/1904/2022"
    mutation = "223V"
    min_tips = 15000
  } else if (subtype == "h1") {
    trim = F
    strain = NA_character_
    mutation = NA_character_
    min_tips = 0
  } else if (subtype == "bvic") {
    trim = T
    strain = "OMAN/2651/2019"
    mutation = "127T"
    min_tips = 15000
  } else if (subtype == "byam") {
    trim = F
    strain = NA_character_
    mutation = NA_character_
    min_tips = 0
  } else {
    stop("Invalid subtype ", subtype)
  }

  # read input files ------------------------------------------------------------

  phy = ape::read.tree(file = inphy)
  fa = seqUtils::fast_fasta(infasta)

  # trim ------------------------------------------------------------
  if (!trim) {
    trimmed_phy = phy
    trimmed_seqs = fa
  } else {
    tree_and_sequences = list(
      tree = phy,
      sequences = tibble::tibble(
        Isolate_unique_identifier = names(fa),
        dna_sequence = seqUtils::clean_sequences(unname(fa), type = "nt")
      )
    )

    usher_tree_and_sequences = convergence::addASRusher(
      tree_and_sequences,
      nuc_ref = fa[1],
      aa_ref = as.character(Biostrings::translate(Biostrings::DNAString(fa[1])))
    )

    usher_tree_and_sequences$tree_tibble$nd = c(
      rep(1, ape::Ntip(usher_tree_and_sequences$tree)),
      castor::count_tips_per_node(usher_tree_and_sequences$tree)
    )

    strain_fullname = usher_tree_and_sequences$tree$tip.label[
      stringr::str_detect(
        usher_tree_and_sequences$tree$tip.label,
        stringr::fixed(strain)
      )
    ]

    parents = rev(convergence:::getParents(
      usher_tree_and_sequences$tree_tibble,
      which(usher_tree_and_sequences$tree_tibble$label == strain_fullname)
    ))

    parents_mutations = usher_tree_and_sequences$tree_tibble$aa_mutations[
      parents
    ]

    ancestor = parents[
      purrr::map_lgl(
        parents_mutations,
        ~ any(stringr::str_detect(.x, stringr::fixed(mutation)))
      )
    ]

    trimmed_phy = castor::get_subtree_at_node(
      usher_tree_and_sequences$tree,
      ancestor - ape::Ntip(usher_tree_and_sequences$tree)
    )$subtree

    if (ape::Ntip(trimmed_phy) < min_tips) {
      stop(
        "Fewer than `min_tips (= ",
        min_tips,
        ") tips in trimmed phylogeny. Perhaps the wrogn branch was found?"
      )
    }

    trimmed_phy = ape::ladderize(trimmed_phy, right = F)
    trimmed_phy$edge.length = trimmed_phy$edge.length / mean(nchar(fa)) # per nt
    trimmed_seqs = fa[names(fa) %in% trimmed_phy$tip.label]
  }

  message("Writing tree")
  ape::write.tree(trimmed_phy, file = outphy)

  message("Writing ", outphy.txt)
  writeLines(outphy, outphy.txt)

  message("Writing ", outfasta)
  seqUtils::write_fast_fasta(
    seqs = unname(trimmed_seqs),
    names = names(trimmed_seqs),
    path = outfasta
  )
  return(0)
}

main()
