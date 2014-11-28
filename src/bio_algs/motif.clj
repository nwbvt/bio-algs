(ns bio-algs.motif
  (:require [bio-algs.ori-rep :refer [varients hamming-dist]]
            [clojure.set :refer :all]))

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

(defn log
  [n base]
  (/ (Math/log n) (Math/log base)))

(defn entropy
  [col]
  (* -1 (apply +
               (for [aa [\A \T \C \G]
                     :let [t (count (filter #(= aa %) col))
                           p (/ t (count col))]]
                 (if (zero? t) 0 (* p (log p 2)))))))

(defn motif-entropy
  "finds the entropy of the given motif"
  [& motifs]
  (let [len (count (first motifs))]
    (apply +
           (for [i (range len) :let [col (map #(nth % i) motifs)]]
             (entropy col)))))

(defn dist
  "Finds the distance between a pattern and one or more strings of dna"
  [pattern dna]
  (if (string? dna)
    (apply min (map (partial hamming-dist pattern) (all-kmers (count pattern) dna)))
    (apply + (map (partial dist pattern) dna))))           ;Sum the dists from each dna string 

(defn median-string
  "Finds the median kmer in the given set of dna
   A median kmer is defined as the one that minimizes the distance between it an all dna strings"
  [k dna]
  (let [kmers (apply union (map (partial all-kmers k) dna))]
    (println kmers)
    (println (map #(dist % dna) kmers))
    (apply min-key #(dist % dna) kmers)))
