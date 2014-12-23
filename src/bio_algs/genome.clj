(ns bio-algs.genome
  (:require [bio-algs.motif :refer [all-kmers]]
            [bio-algs.core :refer [read-file]]
            [clojure.string :refer [join split trim]]))

(defn string-spelled
  "given a sequence of kmers such that each kmer's last k-1 letters
  match the next kmer's first k-1 letter, reconstruct the genome that produced them"
  [seqs]
  (str (apply str (map first seqs)) (apply str (drop 1 (last seqs)))))

(defn overlap-graph
  "returns the overlap graph for a collection of patterns"
  [patterns]
  (for [pre patterns post patterns
        :when (= (drop 1 pre)
                 (take (dec (count post)) post))]
    [pre post]))

(defn format-graph
  "formats the adjacency list in the form of node1 -> node2"
  [adj-list]
  (map #(let [[a b] %] (str a " -> " (if (seq? b) (join "," b) b))) adj-list))

(defn find-indices
  "Returns the indicies in the text where the pattern occurs"
  [text pattern]
  (loop [all [] n (.indexOf text pattern)]
    (if (= -1 n) all
      (recur (conj all n) (.indexOf text pattern (inc n))))))

(defn deBruijn
  "Returns a De Bruijn graph for the given text"
  [k text]
  (let [j (dec k)
        k-1ers (all-kmers j text)]
    (for [pat k-1ers
          :let [indicies (find-indices text pat)
                map-to (for [index indicies :when (< index (- (count text) j))]
                         (.substring text (inc index) (+ (inc index) j)))]
          :when (not-empty map-to)]
      [pat (sort map-to)])))

(defn make-nodes
  "Makes nodes for a deBruijn graph from the edge"
  [edge]
  (if (or (vector? edge) (seq? edge))
    (let [both (map make-nodes edge)] [(map first both) (map second both)])
    [(apply str (take (dec (count edge)) edge)) (apply str (drop 1 edge))]))

(defn deBruijn-kmers
  "Returns a De Bruijn graph from a list of kmers"
  [kmers]
  (apply hash-map (apply concat
   (let [nodes (map make-nodes kmers)
        node-map (group-by first nodes)]
     (for [[pre suf] node-map]
       [pre (map second suf)])))))

(defn read-adj-list
  "Reads in an adjancency list and turns it into a more usable format"
  [adj-list]
  (apply hash-map
         (mapcat #(let [[l r] (split % #"->")]
                    [(trim l) (split (trim r) #",")]) adj-list)))

(defn format-path
  "Formats the path in the desired format"
  [path]
  (join "->" path))

(defn find-cycle
  "Finds a cycle in the graph returns both the cycle and unused edges"
  [graph start]
  (loop [visited [], curr start, nodes graph]
    (let [next-node (nodes curr) new-path (conj visited curr)]
      (if (empty? next-node) [new-path nodes]
        (recur new-path (first next-node) (assoc nodes curr (rest next-node)))))))

(defn find-start
  [graph path]
  (if (empty? path) (first (keys graph))
    (first (for [node path :when (not-empty (graph node))] node))))

(defn insert-in
  [to-insert path]
  (if (empty? path) to-insert
    (let [[pre suf] (split-with #(not= (first to-insert) %) path)]
      (concat pre to-insert (rest suf)))))

(defn euler-cycle
  "Finds an Eulerian cycle given a directed graph"
  [graph & [start]]
  (loop [path [] unused graph max 10 start start]
    (if (some not-empty (vals unused))
      (let [[new-path new-unused] (find-cycle unused (or start (find-start unused path)))]
        (recur (insert-in new-path path) new-unused (dec max) nil))
      path)))

(defn find-unbalanced
  [path]
  (let [freqs-out (fn [k] (count (path k)))
        freqs-in (frequencies (apply concat (vals path)))
        all (set (concat (keys freqs-in) (keys path)))]
    (for [node all
          :let [ins (get-in freqs-in [node] 0)
                outs (freqs-out node)]
          :when (not= ins outs)]
      [(- ins outs) node])))

(defn euler-path
  "Finds an Eulerian path given a directed graph"
  [graph]
  (let [[[_ start] [_ end]] (-> graph find-unbalanced sort)
        g-with-cycle (assoc graph end (conj (get-in graph [end] []) start))
        cyc (euler-cycle g-with-cycle start)]
    (take (dec (count cyc)) cyc)) )

(defn string-recon
  [segments]
  (string-spelled (euler-path (deBruijn-kmers segments))))

(defn cycle-recon
  [segments]
  (string-spelled (euler-cycle (deBruijn-kmers segments))))

(defn bin-string
  [k n]
  (clojure.string/replace (format (str "%" k "s") (Integer/toBinaryString n)) " " "0"))

(defn universal-string
  [k]
  (let [patterns (map (partial bin-string k) (range (Math/pow 2 k)))
        cyc (cycle-recon patterns)
        len (count cyc)]
    ;drop the last kmer
    (apply str (take (- len (dec k)) cyc))))

(defn gapped-string-spelled
  "Spells a string from read pairs"
  [k d seqs]
  (let [f-pats (map first seqs)
        s-pats (map second seqs)
        f-string (string-spelled f-pats)
        s-string (string-spelled s-pats)]
    (apply str (concat (take (+ k d) f-string) s-string))))

(defn recon-from-read-pairs
  "reconstruct a string from read pairs"
  [k d read-pairs]
  (gapped-string-spelled k d (euler-path (deBruijn-kmers read-pairs))))

(defn read-read-pairs
  "Reads in the actual read-pairs from a input file"
  [file]
  (let [raw (read-file file)]
    (map #(split % #"\|") raw)))
