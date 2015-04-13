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

(defn common-start
  "returns the longest common start between two strings"
  [str1 str2]
  (count (take-while #(= (nth str1 %) (nth str2 %)) (range (min (count str1) (count str2))))))

(defn- add-to-suffix-tree
  [full-text old-tree text]
  (loop [tree old-tree parent 0 l text]
    (if (empty? l) tree
      (let [[_ children] (nth tree parent)
            c (first l)
            matching (first (filter #(= c (first (text-at-node (nth tree %) full-text))) children))]
        (if matching
          (let [node-text (text-at-node (nth tree matching) full-text)
                common (common-start l node-text)
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

(defn match-2?
  "Returns whether or not the pattern matches something within the text using a suffix tree"
  [text pattern]
  (let [tree (suffix-tree text)]
    (loop [index 0 l pattern]
      (let [[_ children] (nth tree index)
            match (first (filter #(= (first l) (nth text (first (first (nth tree %))))) children))]
        (if match
          (let [node (nth tree match)
                node-text (text-at-node node text)
                common-count (common-start l node-text)]
            (cond
              (= common-count (count l)) true
              (= common-count (count node-text)) (recur match (drop common-count l))
              :default false))
          false)))))

(defn labels
  "Return the labels of the suffix tree"
  [text tree]
  (for [node tree :when (first node)]
    (text-at-node node text)))

(defn repeats
  "returns the repeats using a suffix tree"
  ([text] (repeats (suffix-tree text) text 0))
  ([tree text index]
   (let [node (nth tree index)]
     (if (empty? (second node))
       []
       (let [repeats-from-children (mapcat (partial repeats tree text) (second node))
             all (conj repeats-from-children "")]
         (map (partial str (if (first node) (text-at-node node text) "")) all))))))

(defn longest-repeat
  "Returns the longest repeat in a string using a suffix tree"
  [text]
  (let [all-repeats (repeats text)] (apply max-key count all-repeats)))

(defn matches-start
  "returns the longest prefix in text that is contained in the suffix tree"
  [suffix-tree full-text to-match]
  (loop [index 0 l to-match prefix ""]
    (let [node (nth suffix-tree index)
          common-node (first (filter #(= (first l) (first (text-at-node (nth suffix-tree %) full-text))) (second node)))]
      (if common-node
        (let [node-text (text-at-node (nth suffix-tree common-node) full-text)
              common-count (common-start l node-text)
              total-common (str prefix (apply str (take common-count l)))]
          (if (= (count node-text) common-count)
            (recur common-node (drop common-count l) total-common)
            total-common))
        prefix))))

(defn longest-common
  "Returns the longest common substring"
  [str1 str2]
  (let [tree (suffix-tree str1)
        match-list (map (partial matches-start tree str1) (take-while identity (iterate next str2)))]
    (apply max-key count match-list)))

(defn- subseqs-of-size
  "Returns all of the subsequences of the text of a particular size"
  [text n]
  (for [start (range (- (count text) n)) :let [end (+ start n)]]
    (subs text start end)))

(defn shortest-unique
  "Returns the shortest non-common substrint"
  [str1 str2]
  (first 
    (filter #(not (match-2? str2 %))
            (mapcat (partial subseqs-of-size str1) (range 1 (count str1))))))

(defn suffix-array
  [text]
  (let [suffixes (for [i (range (count text))] (apply str (drop i text)))]
    (sort-by (partial nth suffixes) (range (count text)))))

(defn bw-transform
  "Transforms into the burrows wheeler transformation of the text"
  [text]
  (let [double-text (concat text text)
        len (count text)
        cycles (map (partial take len)
                    (take len (iterate next double-text)))]
    (apply str (map last (sort-by #(apply str %) cycles)))))

(defn- make-count-list
  [s]
  (loop [counts {} cl [] r s]
    (if (empty? r) cl
      (let [i (first r) c (or (counts i) 0)]
        (recur (assoc counts i (inc c)) (conj cl [i c]) (rest r))))))

(defn bw-recon
  "Transfroms from the burrows wheeler transofrmation back to the text"
  [bw]
  (let [cl (make-count-list bw) 
        sorted (sort cl)
        end (first sorted)]
    (loop [text [] c end]
      (let [i (.indexOf cl c)
            nc (nth sorted i)
            new-text (conj text (first nc))]
        (if (= nc end) (apply str new-text)
          (recur new-text nc))))))

(defn b-search
  "binary search of a sorted vector for x, returning its index"
  [v x]
  (loop [start 0 end (count v)]
    (if (= start end) nil
      (let [mid (int (+ start (/ (- end start) 2)))
            f (nth v mid)
            r (compare f x)]
        (cond
          (neg? r) (recur (inc mid) end)
          (pos? r) (recur start mid)
          :else mid)))))

(defn last-to-first
  "given a burrows wheeler transform, return a vector transforming the index in the
  bw text to its index in the sorted text"
  [bw]
  (let [cl (make-count-list bw)
        sorted (vec (sort cl))]
    (vec (map (partial b-search sorted) cl))))

(defn first-occurences
  "returns a map of the first occurance of a symbol in the sorted version of the given string"
  [input]
  (loop [i 0 res {} s (sort input) prev nil]
    (if (empty? s) res
      (let [c (first s)]
        (if (= prev c) (recur (inc i) res (rest s) prev)
          (recur (inc i) (assoc res c i) (rest s) c))))))

(defn- counts
  [bw]
  (iterate (fn [m]
             (let [i (:i m)
                   c (nth bw i)]
               (assoc m c (inc (m c)) :i (inc i))))
           (zipmap (conj (set bw) :i) (repeat 0))))

(defn count-matrix
  "returns the count matrix for the bw"
  [bw]
  (vec (take (inc (count bw)) (counts bw))))

(defn checkpoint-array
  [bw c]
  (take (/ (inc (count bw)) c) (filter #(zero? (mod (:i %) c)) (counts bw))))

(defn count-from-cp-array
  "returns the count of sym before position i using the bw and checkpoint array"
  [bw cp c sym i]
  (let [cp-i (/ i c)]
   (+ ((nth cp cp-i) sym)
     (count (filter (partial = sym) (subs bw cp-i (+ cp-i (mod i c)))))))) 

(defn- bw-matcher
  "Helper function for bw matching"
  [bw cp-array c pattern]
  (let [fo (first-occurences bw)]
    (loop [start 0 end (dec (count bw)) p (reverse pattern)]
      (if (empty? p) [end start]
        (let [sym (first p)]
          (let [new-start (+ (fo sym) (count-from-cp-array bw cp-array c sym start))
                new-end (dec (+ (fo sym) (count-from-cp-array bw cp-array c sym (inc end))))]
            (recur new-start new-end (rest p))))))))

(defn bw-match-count
  "return the number of times the pattern appears in the burrows wheeler transform"
  ([bw pattern] (bw-match-count bw (count-matrix bw) pattern))
  ([bw counts pattern] (bw-match-count bw counts 1 pattern))
  ([bw cp-array c pattern]
    (inc (apply - (bw-matcher bw cp-array c pattern)))))

(defn bw-match
  "Find the locations of matches to the given patterns"
  [text patterns]
  (let [bw (bw-transform text) 
        suffixes (suffix-array text)]
    
    ))
