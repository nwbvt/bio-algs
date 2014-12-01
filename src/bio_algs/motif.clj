(ns bio-algs.motif
  (:require [bio-algs.ori-rep :refer [varients hamming-dist]]
            [clojure.set :refer :all]
            [clojure.string :refer [split trim]]))

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
    (apply min-key #(dist % dna) kmers)))

(defn score-kmer
  "Scores a kmer against a given profile"
  [kmer profile]
  (apply * 
         (for [i (range (count kmer))
               :let [a (nth kmer i)]]
           (nth (profile a) i))))

(defn most-probable-kmer
  "Finds the most probably kmer in the given dna fitting the given pattern"
  [k dna profile]
  (apply max-key #(score-kmer % profile) (all-kmers k dna)))

(defn make-profile-from-string
  "Makes a profile from the input string format"
  [input]
  (let [strrows (split input #"\n")
        rows (map (fn [r] (map #(Double/parseDouble %) (split (trim r) #" "))) strrows)]
    (zipmap [\A \C \G \T] rows)))

(defn score
  "Score the motifs given"
  [motifs]
  (let [n (count (first motifs))]
   (apply +
         (for [i (range n)
               :let [xs (map (partial #(nth % i)) motifs)
                     freqs (frequencies xs) 
                     consensus (apply max-key freqs (keys freqs))]]
           (count (filter (partial not= consensus) xs))))))

(defn gen-profile
  "Generates a profile from a set of motifs"
  [motifs smooth?]
  (let [acids [\A \C \G \T]
        n (count (first motifs))
        t (float (count motifs))]
    (zipmap acids (for [acid acids]
                    (for [i (range n)
                          :let [total (+ t (if smooth? 4 0))
                                found (+ (count (filter (partial = acid) (map #(nth % i) motifs))) 
                                         (if smooth? 1 0))]]
                      (/ found total))))))

(defn make-motifs
  "the make motif step for greedy motif search"
  [base-motif all-strands k]
  (loop [motifs [base-motif]
         strands all-strands]
    (let [profile (gen-profile motifs true)]
      (if (empty? strands) motifs
        (recur (conj motifs (most-probable-kmer k (first strands) profile)) (rest strands))))))

(defn greedy-motif-search
  "Perform a greedy motif search that returns the best kmer motifs in the given dna strands"
  [k strands]
  (loop [best-motifs (map (partial take k) strands)
         best-score (score best-motifs)
         kmers (all-kmers k (first strands))]
    (if (empty? kmers) best-motifs
      (let [motifs (make-motifs (first kmers) (rest strands) k)
            new-score (score motifs)]
        (if (>= new-score best-score)
          (recur best-motifs best-score (rest kmers))
          (recur motifs new-score (rest kmers)))))))
