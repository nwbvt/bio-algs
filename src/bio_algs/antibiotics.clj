(ns bio-algs.antibiotics
  (:require [clojure.string :refer [split]]
            [bio-algs.ori-rep :refer [reverse-comp]]))

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
  (apply + (map #(mass-table (str %)) peptide)))

(defn theoretical-spectrum
  "Finds the theoretical spectrum of a peptide"
  [peptide]
  (let [pep-len (count peptide)
        pep-cycle (cycle peptide)]
    (sort (conj (for [i (range pep-len) j (range 1 pep-len)
                      :let [sub-pep (->> pep-cycle (drop i) (take j))]]
                  (weight sub-pep)) 0 (weight peptide)))))
