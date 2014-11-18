(ns bio-algs.motif
  (:require [bio-algs.ori-rep :refer [varients]]))

(defn matches?
  "Returns true iff the dna contains the pattern"
  [dna pattern]
  ((comp not neg?) (.indexOf dna pattern)))

(defn matches-with-d-mismatches?
  "returns true iff the dna matches with at most d mismatches"
  [pattern d dna]
  (true? (some (partial matches? dna) (varients pattern d false))))

(defn all-kmers
  "returns all kmers appearing in a given dna sequence"
  [k dna]
  (set (for [r (iterate rest dna) :while (>= (count r) k)]
         (apply str (take k r)))))

(defn motif-enum
  "finds (k,d) motifs in the supplied dna sequences"
  [dna k d]
  (let [dna-source (first dna)
        dna-to-compare (rest dna)]
    (set 
      (for [pattern (all-kmers k (first dna))
            pattern' (varients pattern d false)
            :when (every? (partial matches-with-d-mismatches? pattern' d) dna-to-compare)]
        pattern'))))
