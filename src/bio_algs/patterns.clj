(ns bio-algs.patterns)

(defn- add-trie
  "Adds a pattern to an existing trie"
  [old-trie pattern]
  (loop [trie old-trie parent 0 pat pattern]
    (if (empty? pat) trie
      (let [[_ children] (nth trie parent)
            c (first pat)
            remaining (rest pat)
            matching (first (filter #(= c (first (nth trie %))) children))]
        (if matching (recur trie matching remaining)
          (let [[parent-val siblings] (nth trie parent)
                new-index (count trie)
                new-parent [parent-val (conj siblings new-index)]
                new-entry [c []]]
           (recur (conj (assoc trie parent new-parent) new-entry) new-index remaining)))))))

(defn create-trie
  "Creates a trie from a set of patterns"
  [patterns]
  (reduce add-trie [[nil []]] patterns))

(defn display-trie
  "Displays the trie in the way the course grader wants it"
  [trie]
  (loop [index 0 nodes trie out []]
    (if (empty? nodes) out
      (recur (inc index) (rest nodes)
             (doall (concat out (for [child (second (first nodes))]
                                  (format "%d->%d:%s" index child (first (nth trie child))))))))))

(defn match?
  "Returns whether or not the current text matches a pattern in the trie"
  [trie text]
  (loop [index 0 t text]
    (let [[_ trie-children] (nth trie index)
          c (first t)]
      (if (empty? trie-children) true
        (let [child-vals (map #(first (nth trie %)) trie-children)
              match (.indexOf child-vals c)]
          (if (neg? match) false
            (recur (nth trie-children match) (rest t))) )))))

(defn matches
  "Returns the indexes for which a pattern in the trie matches the text"
  [trie text]
  (filter #(match? trie (drop % text)) (range (count text))))
