(ns bio-algs.antibiotics
  (:require [clojure.string :refer [split join]]
            [bio-algs.ori-rep :refer [reverse-comp]]
            [clojure.math.numeric-tower :refer :all]))

(def rna-codon-table
  (let [text (slurp "resources/RNA_codon_table_1.txt")
        key-vals (map #(split % #" ") (split text #"\n"))]
    (zipmap (map first key-vals) (map second key-vals))))

(def mass-table
  (let [text (slurp "resources/integer_mass_table.txt")
        key-vals (map #(split % #" ") (split text #"\n"))]
    (zipmap (map first key-vals) (map #(Integer/parseInt (second %)) key-vals))))

(defn translate
  "Translate a dna string into the cooresponding amino acid string"
  [dna]
  (let [codons (partition 3 dna)]
    (apply str (for [codon codons :let [codon-string (apply str (replace {\T \U} codon))
                                        amino-acid (rna-codon-table codon-string)]
                     :while amino-acid]
                 amino-acid))))

(defn count-options
  "Counts the number of dna strings that will translate to a given amino acid sequence"
  [& aaseq]
  (apply * (map (frequencies (vals rna-codon-table)) aaseq)))

(defn find-encoding
  "Finds places in the input dna strand that encode to the given amino acid string"
  [dna peptide]
  (let [dna-len (* 3 (count peptide))]
    (loop [dna dna results []]
      (if (< (count dna) dna-len) results
        (let [dna-part (apply str (take dna-len dna))
              forward-pep (translate dna-part)
              reverse-pep (translate (reverse-comp dna-part))
              matches (or (= peptide forward-pep)
                          (= peptide reverse-pep))]
          (recur (rest dna) (if matches (conj results dna-part) results)))))))

(defn weight
  "find the weight of a peptide"
  [peptide]
  (apply + (map #(if (number? %) %
                   (mass-table (str %))) peptide)))

(defn theoretical-spectrum
  "Finds the theoretical spectrum of a peptide"
  [peptide]
  (let [pep-len (count peptide)
        pep-cycle (cycle peptide)]
    (sort (conj (for [i (range pep-len) j (range 1 pep-len)
                      :let [sub-pep (->> pep-cycle (drop i) (take j))]]
                  (weight sub-pep)) 0 (weight peptide)))))

(defn sub-linear
  "counts the number of linear subpeptides for a peptide of length n
   which is basically 1 (for the empty one) + the sum of integers 1-n"
  [n]
  (inc (* n (ceil (/ n 2)))))

(defn bb-seq
  "Uses branch and bounds to find peptide sequences consistent with a spectrum"
  [weights]
  (let [aa-weights (filter (set weights) (set (vals mass-table)))
        weight-set (set weights) ]
    (loop [peptides [[[], 0]] valid []]  ;first value is the sequence (in weights), second the total
      (let [branched (for [[pep-s pep-w] peptides, aa-w aa-weights]
                       [(conj pep-s aa-w) (+ pep-w aa-w) ])
            bounded (filter #(weight-set (second %)) branched)
            valid (concat valid (filter #(= weight-set (set (theoretical-spectrum (first %)))) peptides))]
        (if (empty? bounded) valid
          (recur bounded valid))))))

(defn format-weight-seqs
  "Output the a sequence of weights as a string"
  [weight-seq]
  (join "-" weight-seq))

(defn format-bb-seq
  "Give the bb-seq in the format expected"
  [weights]
  (map format-weight-seqs (map first (bb-seq weights))))


(defn score
  "Scores a spectrum against the theoretical spectrum of a given peptide"
  [peptide spectrum]
  (loop [t-spec (theoretical-spectrum peptide)
         e-spec (vec spectrum)
         c-score 0]
    (if (empty? t-spec) c-score
      (let [w (first t-spec) i (.indexOf e-spec w)]
        (if (= -1 i) (recur (rest t-spec) e-spec c-score)
          (recur (rest t-spec) (assoc e-spec i nil) (inc c-score)))))))
