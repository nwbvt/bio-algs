(ns bio-algs.genome
  (:require [bio-algs.motif :refer [all-kmers]]
            [clojure.string :refer [join]]
            ))

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
  [(apply str (take (dec (count edge)) edge)) (apply str (drop 1 edge))])

(defn deBruijn-kmers
  "Returns a De Bruijn graph from a list of kmers"
  [kmers]
  (sort
   (let [nodes (map make-nodes kmers)
        node-map (group-by first nodes)]
    (for [[pre suf] node-map]
      [pre (sort (map second suf))]))))
