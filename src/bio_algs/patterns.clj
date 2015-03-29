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

(defn- text-at-node
  [[[start len] _] text]
  (subs text start (+ start len)))

(defn- add-to-suffix-tree
  [full-text old-tree text]
  (loop [tree old-tree parent 0 l text]
    (if (empty? l) tree
      (let [[_ children] (nth tree parent)
            c (first l)
            matching (first (filter #(= c (first (text-at-node (nth tree %) full-text))) children))]
        (if matching
          (let [node-text (text-at-node (nth tree matching) full-text)
                common (count (take-while #(= (nth l %) (nth node-text %)) (range (min (count l) (count node-text)))))
                new-entry [[(+ (- (count full-text) (count l)) common) (- (count l) common)] []]]
           (if (= common (count node-text))
            (recur tree matching (drop common l))
            (let [common (count (take-while #(= (nth l %) (nth node-text %)) (range (count l))))
                  [[start old-end] siblings] (nth tree matching)
                  new-matching [[start common]
                                [(count tree) (inc (count tree))]]
                  extended [[(+ start common) (- old-end common)] siblings]]
              (conj (assoc tree matching new-matching) extended new-entry))))
          (let [[parent-val siblings] (nth tree parent)
                new-parent [parent-val (conj siblings (count tree))]
                new-entry [[(- (count full-text) (count l)) (count l)] []]]
            (conj (assoc tree parent new-parent) new-entry)))))))

(defn suffix-tree
  "Creates a suffix tree for the given text"
  [text]
  (reduce (partial add-to-suffix-tree text)
          [[nil []]]
          (take-while identity (iterate next text))))

(defn labels
  "Return the labels of the suffix tree"
  [text tree]
  (for [node tree :when (first node)]
    (text-at-node node text)))
