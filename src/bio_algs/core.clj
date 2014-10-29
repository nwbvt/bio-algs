(ns bio-algs.core)

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

(defn most-common-kmer
  "find the most common kmer in the text"
  [text k]
  (loop [left text, counts {}]
    (if (>= (count left) k)
      (let [part (take k left)
            new-count (inc (get-in counts [part] 0))]
        (recur (rest left) (assoc counts part new-count)))
      (let [max-count (apply max (vals counts))]
        (for [[p c] counts :when (= c max-count)] (apply str p))))))

(defn reverse-comp
  "find the reverse complement of the given dna string"
  [dna]
  (let [comps {\A \T, \T \A, \C \G, \G \C}]
    (apply str (map comps (reverse dna)))))

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
  (loop [left text kmers #{} counts {} drop-end (concat (repeat L nil) text)]
    (if (< (count left) k) (set (map #(apply str %) kmers))
      (let [add-kmer (take k left)
            rem-kmer (take k drop-end)
            new-counts (assoc 
                         (assoc counts 
                                rem-kmer
                                (dec (get-in counts [rem-kmer] 0)))
                         add-kmer 
                         (inc (get-in counts [add-kmer] 0)))
            meets-crit (= (new-counts add-kmer) t)]
        (recur (rest left) (if meets-crit (conj kmers add-kmer) kmers) new-counts (rest drop-end))))))
