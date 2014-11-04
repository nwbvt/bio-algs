(ns bio-algs.antibiotics
  (:require [clojure.string :refer [split]]))

(def rna-codon-table
  (let [text (slurp "resources/RNA_codon_table_1.txt")
        key-vals (map #(split % #" ") (split text #"\n"))]
    (zipmap (map first key-vals) (map second key-vals))))

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
