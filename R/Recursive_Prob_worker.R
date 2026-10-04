
#processPed()
#configProb_recursive_oneFounder()
#computeSharingProb_recursive()

# 2. Recursive familial configuration probability workers
###############################################################################

processPed <- function(ped, carriers = NULL) {
  if ("Gen" %in% colnames(ped) & any(is.na(ped$Gen))) {
    ped <- ped[, -which(colnames(ped) == "Gen")]
  }

  ped <- RVMethods:::affected_onlyPed(ped)
  ped <- ped[order(ped$ID), ]

  pdgree <- SimRVPedigree::ped2pedigree(ped)
  origID <- pdgree$id
  pdgree$id <- seq_along(pdgree$id)

  parents <- sapply(pdgree$id, function(i) c(pdgree$findex[i], pdgree$mindex[i]))
  founders <- which(parents[1, ] == 0)
  finalDescendants <- which(!(pdgree$id %in% c(parents[1, ], parents[2, ])))

  if (inherits(pdgree$affected, "matrix")) {
    pdgree$affected <- pdgree$affected[, 1]
  }

  affected <- which(pdgree$affected == 1)

  if (!is.null(carriers)) {
    carriers <- which(origID %in% carriers)
  } else {
    carriers <- affected
  }

  if (length(affected) < 2) stop("need at least 2 affected subjects")
  if (sum(!carriers %in% affected) > 0) stop("carriers must be a subset of affected")

  list(
    origID = origID,
    ped = pdgree,
    parents = parents,
    founders = founders,
    affected = affected,
    size = length(pdgree$id),
    carriers = carriers,
    id = pdgree$id,
    finalDescendants = finalDescendants
  )
}

# configProb_recursive_oneFounder

configProb_recursive_oneFounder <- function(procPed, founder, tau_vec) {
  n <- procPed$size
  parents <- procPed$parents
  affected <- procPed$affected
  carriers <- procPed$carriers

  tau_vec <- as.numeric(tau_vec)
  nt <- length(tau_vec)

  children <- vector("list", n)
  for (j in seq_len(n)) children[[j]] <- integer(0)

  for (child in seq_len(n)) {
    pa <- parents[, child]
    pa <- pa[pa > 0]
    for (p in pa) {
      children[[p]] <- c(children[[p]], child)
    }
  }

  visit <- function(j, postorder) {
    for (k in children[[j]]) {
      postorder <- visit(k, postorder)
    }
    c(postorder, j)
  }

  postorder <- visit(founder, integer(0))

  carrier_not_desc <- carriers[!(carriers %in% postorder)]
  if (length(carrier_not_desc) > 0) {
    return(rep(0, nt))
  }

  p0 <- vector("list", n)
  p1 <- vector("list", n)

  for (j in postorder) {
    if (j %in% affected) {
      if (j %in% carriers) {
        I0 <- 0; I1 <- 1
      } else {
        I0 <- 1; I1 <- 0
      }
    } else {
      I0 <- 1; I1 <- 1
    }

    pj0 <- rep(I0, nt)
    pj1 <- rep(I1, nt)

    for (k in children[[j]]) {
      pj0 <- pj0 * p0[[k]]
      pj1 <- pj1 * ((1 - tau_vec) * p0[[k]] + tau_vec * p1[[k]])
    }

    p0[[j]] <- pj0
    p1[[j]] <- pj1
  }

  p1[[founder]]
}

# computeSharingProb_recursive

computeSharingProb_recursive <- function(procPed, tau_vec, pc) {
  tau_vec <- as.numeric(tau_vec)

  nfounders <- length(procPed$founders)
  out <- rep(0, length(tau_vec))

  for (f in procPed$founders) {
    out <- out + configProb_recursive_oneFounder(
      procPed = procPed,
      founder = f,
      tau_vec = tau_vec
    )
  }

  zeroconfig <- as.numeric(length(procPed$carriers) == 0)

  pc * out + zeroconfig * (1 - nfounders * pc)
}
