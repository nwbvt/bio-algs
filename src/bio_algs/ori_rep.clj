(ns bio-algs.ori-rep
  (:use clojure.java.io clojure.set))

(defn pattern-count
  "count the number of times the pattern occurs in the text"
  [text pattern]
  (let [pat-len (count pattern)
        pat-seq (seq pattern)]
    (loop [pat-count 0
           rem-text text]
      (if (< (count rem-text) pat-len)
        pat-count
        (let [matches (= (take pat-len rem-text) pat-seq)]
         (recur (if matches (inc pat-count) pat-count) (rest rem-text))))))) 

(defn swap
  "returns the map passed in with the given function applied to the given key"
  [m f k & [default]]
  (let [old (get-in m [k] default)]
    (assoc m k (f old))))

(defn reverse-comp
  "find the reverse complement of the given dna string"
  [dna]
  (let [comps {\A \T, \T \A, \C \G, \G \C}]
    (apply str (map comps (reverse dna)))))

(defn- mutate
  "returns the single letter mutations of a given dna sequence"
  [dna]
  (let [dna-v (apply vector dna)]
    (map #(apply str %)
      (set (apply concat
        (for [i (range (count dna))]
          (map #(assoc dna-v i %) [\A \T \C \G])))))))

(defn varients
  "returns the varients of the DNA sequence that are witing the given hamming distance"
  [dna d include-reverse?]
  (let [mutations
    (loop [vars #{dna} left d]
      (if (zero? left) vars
        (recur (set (apply union (conj (map mutate vars) vars))) (dec left))))]
    (if include-reverse? (concat mutations (map reverse-comp mutations))
      mutations)))

(defn most-common-kmer
  "find the most common kmer in the text, optionally including a minimum hamming distance for approx matches"
  ([text k] (most-common-kmer text k 0 false))
  ([text k d] (most-common-kmer text k d false))
  ([text k d include-reverse?]
   (let [varient-parts (apply concat (for [left (iterate rest text) :let [part (apply str (take k left))] 
                                           :while (>= (count left) k)]
                                       (varients part d include-reverse?)))
         counts (frequencies varient-parts)
         max-count (apply max (vals counts))]
     (for [[p c] counts :when (= c max-count)] p))))

(defn pat-match
  "find an exact pattern match"
  [pattern text]
  (let [len (count pattern)]
    (loop [left text matches [] loc 0]
      (if (< (count left) len)
        matches
        (recur (rest left)
               (if (= (seq pattern) (take len left))
                 (conj matches loc)
                 matches)
               (inc loc))))))

(defn find-clumps
  "Find all distinct k-mers forming (L,t)-clumps in the text
   k is the length of the k-mer, L is the length of text they clump exists in,
   and t is the amount of times the kmer needs to occur"
  [text k L t]
  (loop [left text drop-end (concat (repeat (inc (- L k)) nil) text) kmers #{} counts {}]
    (if (< (count left) k) (set (map #(apply str %) kmers))
      (let [add-kmer (take k left)
            rem-kmer (take k drop-end)
            new-counts (-> counts (swap dec rem-kmer 0)  ;Remove the kmer moving out of the window  
                                  (swap inc add-kmer 0)) ;Add the kmer moving in the window   
            meets-crit (= (new-counts add-kmer) t)] ;If the count of the kmer just added is at our threshold 
        (recur (rest left) (rest drop-end) 
               (if meets-crit (conj kmers add-kmer) kmers)
               new-counts)))))

(let [skew-map {\A 0, \T 0, \C -1, \G 1}]
  (defn skew
    "Returns the skew (difference between the count of gs and cs) of the dna sequence"
    [dna]
    (loop [rem-dna dna, skews [0]]
      (if (empty? rem-dna) skews
        (recur (rest rem-dna) (conj skews (+ (last skews) (skew-map (first rem-dna))))))))

  (defn min-skew
  "Returns the indices where the skew hits a minimum"
  [dna]
  (loop [rem-dna dna
         min-val 0
         min-index []
         cur-val 0
         cur-index 0]
    (if (empty? rem-dna) min-index
      (let [new-val (+ cur-val (skew-map (first rem-dna)))
            min-diff (- min-val cur-val)]
       (recur (rest rem-dna)
              (if (pos? min-diff) cur-val min-val)
              (if (pos? min-diff) [cur-index] (if (zero? min-diff) (conj min-index cur-index) min-index))
              new-val
              (inc cur-index)))))))


(defn hamming-dist
  "returns the hamming distance between two sequences"
  [& seqs]
  (let [len (count seqs)]
    (count (filter #(not (apply = %)) (partition len (apply interleave seqs))))))
