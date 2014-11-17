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
  [aaseq]
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
  ([peptide] (theoretical-spectrum peptide false))
  ([peptide linear?]
  (sort
    (let [pep-len (count peptide)
          premasses (loop [prev 0 masses [0] peps peptide]
                      (if (empty? peps) masses
                        (let [next-weight (+ prev (weight [(first peps)]))]
                          (recur next-weight (conj masses next-weight) (rest peps)))))]
      (concat
        [0]
        (for [i (range pep-len) j (range i pep-len)]
          (- (premasses (inc j)) (premasses i)))
        (if linear? []
          (for [i (range 1 pep-len) j (range i (dec pep-len))]
            (- (premasses pep-len) (- (premasses (inc j)) (premasses i))))))))))

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
  ([peptide spectrum] (score peptide spectrum false))
  ([peptide spectrum linear?]
  (loop [t-spec (theoretical-spectrum peptide linear?)
         e-spec (if (map? spectrum) spectrum ; if it is already in the form of the frequency map
                  (frequencies spectrum))
         c-score 0]
    (if (empty? t-spec) c-score
      (let [w (first t-spec) i (e-spec w)]
        (if (and (not (nil? i)) (not (zero? i))) (recur (rest t-spec) (assoc e-spec w (dec i)) (inc c-score)) 
          (recur (rest t-spec) e-spec c-score)))))))

(defn trim
  "trims the leaderboard"
  [leaderboard n]
  (let [sorted (reverse (sort-by :score leaderboard))
        best (take n sorted)
        remain (drop n sorted)
        worst-allowed (:score (last best))
        tied (take-while #(= worst-allowed (:score %)) remain)]
    (concat best tied)))

(defn lb-sequence
  "Uses the leaderboard algorithm to find the peptide sequence"
  ([n weights] (lb-sequence n weights (vals mass-table)))
  ([n weights aas]
  (let [aa-weights (set aas)
        parent-mass (apply max weights)
        weight-freqs (frequencies weights)
        mk-entry (fn [pep] {:pep pep :mass (apply + pep) :score (score pep weight-freqs true)})]
    (loop [peptides [[]] leader (mk-entry [])]
      (let [branched (for [pep peptides, aa aa-weights :let [new-pep (conj pep aa)]]
                       (mk-entry new-pep))
            potential (filter #(= parent-mass (:mass %)) branched)
            best-potential (first (reverse (sort-by :score (conj potential leader))))
            left (filter #(> parent-mass (:mass %)) branched) ]
        (if (empty? left) (:pep best-potential)
          (let [best-left (trim left n)]
            (recur (map :pep best-left) best-potential))))))))


(defn convolution
  "Computes the convolution of the spectrum"
  [spectrum]
  (let [sorted (vec (sort spectrum))
        spec-count (count spectrum) ]
    (for [i (range spec-count) j (range (inc i) spec-count)
          :let [ith (sorted i) jth (sorted j)] :when (not (= ith jth))]
      (- jth ith))))

(defn best-convolutions
  "Gets the most frequent m convolutions"
  [m spectrum]
  (let [convolutions (convolution spectrum)
        freqs (frequencies convolutions)
        sorted (reverse (sort-by second freqs))]
        (map first (take m sorted))) )

(defn convolution-seq
  "Uses the most frequent convolution masses as the basis for a leaderboard sequence"
  [m n spectrum]
  (let [best (best-convolutions m spectrum)]
    (lb-sequence n spectrum best)))
