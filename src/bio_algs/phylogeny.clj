(ns bio-algs.phylogeny
  (:require [clojure.string :refer [split]]
            [clojure.set :refer [intersection difference]]
            [bio-algs.core :refer [read-file to-ints write-result parse-graph]]
            [bio-algs.ori-rep :refer [hamming-dist]]))

(defn- is-leaf?
  "Returns true iff the given node is a leaf for the given graph"
  [graph node]
  (= 1 (count (graph node))))

(defn- walk-path
 "find the path and distance between two nodes"
 [graph node1 node2]
 (loop [nodes {node1 [[node1] 0.0]}]
   (if (nil? (nodes node2))
     (if (empty? nodes) [[] Double/POSITIVE_INFINITY]
       (recur (apply merge (for [[node [path d]] nodes
                                 :let [children (apply dissoc (graph node) (keys nodes))
                                       children-nodes (keys children)]]
                             (zipmap children-nodes (map #(vector (conj path %) (+ d (children %))) children-nodes))))))
     (nodes node2))))

(defn- distance
  "find the distance between two nodes in a graph"
  [graph node1 node2]
  (second (walk-path graph node1 node2)))

(defn- path
  "Find the path from one node to another"
  [graph node1 node2]
  (first (walk-path graph node1 node2)))

(defn weight-matrix
  "Find the distance matrix of all leaves to the other leaves"
  [graph]
  (let [leaves (sort (filter (partial is-leaf? graph) (keys graph)))]
    (vec (for [leaf1 leaves]
           (vec (for [leaf2 leaves] (distance graph leaf1 leaf2)))))))

(defn limb-length
  "Find the limb length of a node to its parent for graph with the given distance matrix"
  [matrix i]
  (let [n (count matrix)]
    (apply min (for [j (range n) k (range n) :when (and (not= i j) (not= i k))]
                 (/ (- (+ (-> matrix (nth i) (nth j)) (-> matrix (nth i) (nth k))) (-> matrix (nth j) (nth k))) 2)))))

(defn- find-insertion-spot
  "finds the insertion spot for the limb for the last leaf
  It will find two leaves i,j such that Dij = Din + Djn"
  [matrix]
  (let [n (dec (count matrix))]
    (first
      (for [i (range n) j (range n) 
            :when (= (get-in matrix [i j]) (+ (get-in matrix [i n]) (get-in matrix [j n])))]
        [i j]))))

(defn find-spot-in-path
  "finds the spot in the path the given distance away"
  [graph distance path]
  (loop [p path d distance]
    (let [a (first p) b (second p)
          part ((graph a) b)]
      (assert (not (empty? p)))
      (if (>= d part)
        (recur (rest p) (- d part))
        [[a b] d]))))

(defn- adaptive-phylo
  "helper function for adaptive phylogeny (non tail) recursive algorithm"
  [matrix original-n]
  (let [n (dec (count matrix))]
    (if (= 1 n) {0 {1 (get-in matrix [0 1])}
                 1 {0 (get-in matrix [1 0])}} 
      (let [limb (limb-length matrix n)
            nth-row (vec (map #(- % limb) (nth matrix n)))
            bald-partial (assoc matrix n nth-row)
            bald (vec (map-indexed (fn [i row] (assoc row n (nth-row i))) bald-partial))
            [node1 node2] (find-insertion-spot bald)
            new-matrix (subvec (vec (map #(subvec % 0 n) bald)) 0 n)
            t (adaptive-phylo new-matrix original-n)
            path (path t node1, node2)
            [[i j] d-to-i] (find-spot-in-path t (- ((matrix node1) n) limb) path)]
        (if (zero? d-to-i)
          ;just need to append to existing inner node i
          (-> t (assoc-in [i n] limb)
                (assoc-in [n i] limb))
          ;need to create a new innter node d away from i between i and j
          (let [next-node (+ original-n (- (count (keys t)) n))
                old-i (t i) old-j (t j) old-d ((t i) j)
                d-to-j (- old-d d-to-i)
                new-i (assoc (dissoc old-i j) next-node d-to-i)
                new-j (assoc (dissoc old-j i) next-node d-to-j)]
            (assoc t i new-i j new-j next-node {i d-to-i j d-to-j n limb} n {next-node limb}))))))) 

(defn make-tree
  "Make the tree matching the given distance matrix"
  [matrix]
  (adaptive-phylo matrix (count matrix)))

(defn format-graph
  "Formats the graph"
  [graph]
  (mapcat (fn [node] (map (fn [[k v]] (if (nil? v) (format "%s->%s" node k)
                                        (format "%s->%s:%f" node k v))) (graph node))) (keys graph)))

(def cluster-dist
  "Find the distance between two clusters given a distance matrix"
  (memoize
    (fn [d-matrix c1 c2]
      (/ (apply + (for [m1 (:members c1) m2 (:members c2)]
                    ((d-matrix m1) m2)))
         (* (count (:members c1)) (count (:members c2)))))))

(defn- nearest-clusters
  "Find the clusters that are nearest given a distance matrix"
  [clusters d-matrix]
  (apply min-key first
         (for [c1 clusters c2 clusters :when (not= c1 c2)] [(cluster-dist d-matrix c1 c2) c1 c2])))

(defn upgma
  "Creates a ultrametric graph using unweighted pair group method with arithmetic mean"
  [matrix]
  (let [n (count matrix)]
   (loop [clusters (set (map #(hash-map :index % :members (list %)) (range n)))
          next-node n
          ages (vec (repeat n 0))
          t (zipmap (range n) (repeat n {}))]
     (if (= (count clusters) 1) t
       (let [[dist c1 c2] (nearest-clusters clusters matrix)
             age (/ dist 2)
             weight1 (double (- age (ages (:index c1))))
             weight2 (double (- age (ages (:index c2))))]
         (recur (conj (disj clusters c1 c2)
                      {:index next-node :members (concat (:members c1) (:members c2))})
                (inc next-node)
                (conj ages age)
                (-> t
                    (assoc-in [(:index c1) next-node] weight1)
                    (assoc-in [next-node (:index c1)] weight1)
                    (assoc-in [(:index c2) next-node] weight2)
                    (assoc-in [next-node (:index c2)] weight2))))))))

(defn- find-totals
  [d]
  (vec (map #(apply + %) d)))

(defn d*
  "create the D* matrix for neighbor joining algorithm"
  ([matrix] (d* matrix (find-totals matrix)))
  ([matrix totals]
   (let [n (count matrix)
         nodes (range n)]
     (vec (for [n1 nodes]
            (vec (for [n2 nodes]
                   (if (= n1 n2) 0
                     (- (* (- n 2) ((matrix n1) n2)) (totals n1) (totals n2))))))))))

(defn min-val
  "Find the coordinates of the given matrix that is the minimum"
  [matrix]
  (let [nodes (range (count matrix))]
    (apply min-key (fn [[i j]]
                     ((matrix i) j))
           (for [i nodes j nodes :when (not= i j)] [i j]))))

(defn remove-from-vec
  [v & indices]
  (vec (for [i (range (count v)) :when (neg? (.indexOf indices i))]
         (v i))))

(defn- make-new-d
  [d i j]
  (let [n (count d)
        dij ((d i) j)
        new-dists (conj (mapv (fn [k]
                                (/ (- (+ ((d k) i) ((d k) j)) dij) 2))
                              (remove-from-vec (vec (range n)) i j))
                        0)
        pruned (mapv #(remove-from-vec % i j) (remove-from-vec d i j))]
    (conj
      (vec (for [i (range (- n 2))]
             (conj (pruned i) (new-dists i))))
      new-dists)))

(defn nj-tree
  "Create a graph from the given distance matrix using the neighbor joining algorithm"
  ([d] (nj-tree d (vec (range (count d))) (count d)))
  ([d labels] (nj-tree d labels 0))
  ([d labels next-node]
   (let [n (count d)]
     (assert (< 1 n) (str "distance matrix " d " too small"))
     (if (= 2 n)
       (let [[n1 n2] labels
             dist (double ((d 1) 0))]
         {n1 {n2 dist}
          n2 {n1 dist}})
       (let [totals (find-totals d)
             [i j] (min-val (d* d totals))
             t (nj-tree (make-new-d d i j)
                        (conj (remove-from-vec labels i j) next-node)
                        (inc next-node))
             delta (/ (- (totals i) (totals j)) (- n 2))
             dij ((d i) j)
             i-node (labels i)
             j-node (labels j)
             i-limb (double (/ (+ dij delta) 2))
             j-limb (double (/ (- dij delta) 2))]
         (-> t
             (assoc-in [next-node i-node] i-limb)
             (assoc-in [i-node next-node] i-limb)
             (assoc-in [next-node j-node] j-limb)
             (assoc-in [j-node next-node] j-limb)))))))

(defn run-nj-tree!
  [infile]
  (let [d (to-ints (rest (read-file infile)))
        result (nj-tree d)]
    (write-result (format-graph result))))

(defn tree-leaves
  "Returns the leaves of a tree"
  [tree]
  (let [in-nodes (set (keys tree))
        out-nodes (set (mapcat keys (vals tree)))]
    (difference out-nodes in-nodes)))

(defn next-ripe
  "Find the next ripe node"
  [tree tags]
  (let [ripe-nodes (filter (fn [node] (and (not (:tag (tags node)))
                                           (every? #(:tag (tags %)) (keys (tree node)))))
                           (keys tree))]
   (first ripe-nodes)))

(defn- cost-for-k
  [tags k child]
  (let [c-vals (get-in tags [child :c-vals])]
    (apply min (map #(+ (if (= % k) 0 1) (c-vals %))
                    (keys c-vals)))))

(defn- next-cost-map
  [tags children alphabet]
  (zipmap alphabet (map (fn [k] (apply + (map (partial cost-for-k tags k) children))) alphabet)))

(defn- backtrack-tags
  [tree tags node & [parent]]
  (let [c-vals (get-in tags [node :c-vals])
        alphabet (keys c-vals)
        unchanged-cost (c-vals parent)
        min-choice (apply min-key c-vals alphabet)
        min-cost (c-vals min-choice)
        choice (if (or (nil? unchanged-cost) (<= 1 (- unchanged-cost min-cost)))
                 min-choice
                 parent)]
    (if (nil? (tree node)) {node choice}
      (let [children (keys (tree node))
            child-choices (apply merge (map #(backtrack-tags tree tags % choice) children))]
        (assoc child-choices node choice)))))

(defn partial-sp
  "partial method for running small parsimoney with single character leaves
  the utilized character will be the ith character in each leaf
  This will return a mapping from each node to a seq consisting of its string and the distance to it"
  [tree alphabet leaves i]
  (loop [tags (merge
                (zipmap (keys tree) (repeat {:tag false}))
                (zipmap leaves (for [leaf leaves] {:tag true :c-vals {(nth leaf i) 1}})))
         last-node nil]
    (let [ripe (next-ripe tree tags)]
      (if (nil? ripe) (backtrack-tags tree tags last-node)
        (recur (assoc tags ripe {:tag true :c-vals (next-cost-map tags (keys (tree ripe)) alphabet)}) ripe)))))

(defn- combine-tree-slices
  [slices nodes]
  (zipmap nodes (map (fn [node] (apply str (map #(% node) slices))) nodes)))

(defn small-parsimony-mapping
  "Runs the small parsimony algorithm to find ancestor sequences"
  ([tree] (small-parsimony-mapping tree #{\A \T \C \G}))
  ([tree alphabet]
   (let [leaves (tree-leaves tree)
         parts (map (partial partial-sp tree alphabet leaves) (range (count (first leaves))))]
     (combine-tree-slices parts (keys tree)))))

(defn total-cost
  "Finds the total cost of a small parsimony mapping for a given tree"
  [tree mapping]
  (apply + (for [[node children] tree]
             (apply + (for [child (keys children)]
                        (hamming-dist (mapping node) (if (string? child) child (mapping child))))))))

(defn graph-edges-with-mapping
  "Retuns string versions of the graph (with edges made bidirectional) with a given mapping from a unidirectional tree"
  [tree mapping]
  (mapcat (fn [[node children]]
            (apply concat
                   (for [child (keys children)]
                     (let [child-string (if (string? child) child (mapping child))
                           node-string (mapping node)
                           dist (hamming-dist node-string child-string)]
                        [(str node-string "->" child-string ":" dist)
                         (str child-string "->" node-string ":" dist)]))))
          tree))

(defn root-tree
  "Creates a rooted unidirectional tree from an unrooted bidrectional graph"
  ([tree] (root-tree tree (first (filter integer? (keys tree)))))
  ([tree root]
   (let [children (keys (tree root))]
     (if (empty? children) {} ; hit the leaf
       (let [cleansed (reduce (fn [t child]
                                (assoc t child (dissoc (t child) root)))
                              tree children)]
         (apply merge {root (zipmap children (map (tree root) children))} (map (partial root-tree cleansed) children))))))) 

(defn run-small-parsimony!
  [input-file rooted?]
  (let [input (read-file input-file)
        raw-tree (parse-graph (rest input))
        tree (if rooted? raw-tree (root-tree raw-tree))
        mapping (small-parsimony-mapping tree)]
    (concat [(total-cost tree mapping)] (graph-edges-with-mapping tree mapping))))

(defn nearest-neighbors
  [node1 node2 tree]
  (let [node1-children (keys (tree node1))
        node2-children (keys (tree node2))]
    (assert (= 3 (count node1-children) (count node2-children)))
    (assert (.contains node1-children node2))
    (assert (.contains node2-children node1))
    (let [[w x] (remove #(= node2 %) node1-children)
          [y z] (remove #(= node1 %) node2-children)
          [w-children x-children y-children z-children] (map tree [w x y z])
          swap (fn [children in out] (assoc (dissoc children out) in nil))]
      [(assoc tree node1 {node2 nil w nil y nil} node2 {node1 nil x nil z nil} y (swap y-children node1 node2) x (swap x-children node2 node1))
       (assoc tree node1 {node2 nil w nil z nil} node2 {node1 nil x nil y nil} z (swap z-children node1 node2) x (swap x-children node2 node1))])))

(defn run-nearest-neighbors!
  [input-file]
  (let [input (read-file input-file)
        [node1 node2] (map #(Integer/parseInt %) (split (first input) #" "))
        tree (parse-graph (rest input))
        [tree1 tree2] (nearest-neighbors node1 node2 tree)]
    (write-result (concat (format-graph tree1) [""] (format-graph tree2)))))

(defn neighbor-pairs
  "lists all the neighbor pairs that do not include leaves in a rooted tree"
  [tree]
  (let [internal-nodes (keys tree)]
    (for [node1 internal-nodes node2 (keys (tree node1)) :when (.contains internal-nodes node2)]
      [node1 node2])))

(defn large-parsimony
  "Use the nearest neighbors heuristic to solve the large parismony problem given an unrooted graph
   This will return a seq of each iteration"
  [graph]
  (loop [current graph
         rooted (root-tree graph)
         cur-mapping (small-parsimony-mapping rooted)
         cur-score (total-cost rooted cur-mapping)
         prior []]
    (let [neighbors (vec (mapcat (fn [[node1 node2]]
                                   (nearest-neighbors node1 node2 current))
                                 (neighbor-pairs rooted)))
          n (count neighbors)
          rooted-neighbors (vec (map root-tree neighbors))
          mappings (vec (map small-parsimony-mapping rooted-neighbors))
          costs (vec (map #(total-cost (nth rooted-neighbors %) (nth mappings %)) (range n)))
          best-index (apply min-key costs (range n))
          state {:cost cur-score :mapping cur-mapping :graph current}]
      (if (<= cur-score (costs best-index)) (conj prior state)
        (recur (neighbors best-index) (rooted-neighbors best-index) (mappings best-index) (costs best-index) (conj prior state))))))

(defn run-large-parsimony
  [input-file]
  (let [input (rest (read-file input-file))
        results (large-parsimony (parse-graph input))]
    (mapcat #(concat [(:cost %)] (graph-edges-with-mapping (root-tree (:graph %)) (:mapping %)) [""]) results)))
