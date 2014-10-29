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
