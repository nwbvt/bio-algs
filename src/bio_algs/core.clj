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

