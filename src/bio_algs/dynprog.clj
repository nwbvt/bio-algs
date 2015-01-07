(ns bio-algs.dynprog
  (:require [clojure.string :refer [split]])
  (:require [clojure.set :refer :all]))

(defn find-change
  "given a money value and a list of coins, return the minimum number of coins needed to make that amount"
  [value coins]
  (loop [min-nums [0] v 1]
    (let [m (inc (apply min (for [c coins :when (<= c v)]
                              (min-nums (- v c)))))]
      (if (= value v) m
        (recur (conj min-nums m) (inc v))))))

(defn manhattan
  "Finds the longest distance in the manhatten tourist problem"
  [n m down right]
  (let [s (atom {})
        helper (fn [s i j]
                 (let [d-val (if (zero? i) 0 (+ (nth (nth down (dec i)) j) (s [(dec i) j])))
                      r-val (if (zero? j) 0 (+ (nth (nth right i) (dec j)) (s [i, (dec j)])))
                      v (max d-val r-val)]
                  (assoc s [i,j] v)))]
    (doseq [i (range 0 (inc n)) j (range 0 (inc m))]
      (swap! s helper i j))
    (@s [n, m])))

(defn read-in-manhattan
  [file]
  (let [make-ints (fn [strings] (map #(Integer/parseInt %) strings))
        lines (map #(make-ints (split % #" ")) (split (slurp file) #"\n"))
        [n, m] (first lines)
        down (take n (drop 1 lines))
        right (drop (inc n) lines)]
    {:n n :m m :down down :right right}))

(defn find-backtrack
  "backtracks through the graph to find the longest common subseq"
  [s seq1 seq2]
  (loop [i (count seq1) j (count seq2) align1 '() align2 '()]
    (if (and (zero? i) (zero? j)) [(:v (s [(count seq1) (count seq2)])) (apply str align1) (apply str align2)]
      (let [v (cond (zero? i) {:d :right} (zero? j) {:d :down} :else (s [i, j]))]
        (case (:d v)
          :down   (recur (dec i) j (conj align1 (nth seq1 (dec i))) (conj align2 "-"))
          :right  (recur i (dec j) (conj align1 "-") (conj align2 (nth seq2 (dec j))))
          :across (recur (dec i) (dec j) (conj align1 (nth seq1 (dec i))) (conj align2 (nth seq2 (dec j)))))))))

(defn longest-common-subseq
  "find the longest common subsequence of two dna strands"
  ([seq1 seq2] (longest-common-subseq seq1 seq2 0))
  ([seq1 seq2 indel-pen] (longest-common-subseq seq1 seq2 indel-pen (fn [i j] (if (= i j) 1 0))))
  ([seq1 seq2 indel-pen score-fn]
    (let [s (atom {[0 0] {:v 0}})
          helper (fn [s i j]
                   (assoc s [i,j]
                     (let [d-val {:v (if (zero? i) Double/NEGATIVE_INFINITY (- (:v (s [(dec i) j])) indel-pen)) :d :down}
                           r-val {:v (if (zero? j) Double/NEGATIVE_INFINITY (- (:v (s [i (dec j)])) indel-pen)) :d :right}
                           across-val {:v (if (or (zero? i) (zero? j)) Double/NEGATIVE_INFINITY
                                            (+ (score-fn (nth seq1 (dec i)) (nth seq2 (dec j))) 
                                               (:v (s [(dec i) (dec j)]))))
                                       :d :across}]
                       (max-key :v d-val r-val across-val))))]
      (doseq [i (range 0 (inc (count seq1))) j (range 0 (inc (count seq2))) :when (or (pos? i) (pos? j))]
        (swap! s helper i j))
      (find-backtrack @s seq1 seq2))))

(defn top-order
  [graph]
  (let [all-edges (map #(vec [(first %) (second %)]) graph)
        all-nodes (set (apply concat all-edges))]
    (loop [order []
           edges (set all-edges)]
      (let [out-nodes (set (map second edges))
            candidates (difference all-nodes out-nodes (set order))]
        (if (empty? candidates)
          (do (assert (empty? edges) "Graph was not a DAG") order)
          (let [n (first candidates)
                edges-to-remove (filter #(= n (first %)) edges)]
            (recur (conj order n) (difference edges edges-to-remove))))))))

(defn longest-common-path
  "solves a generic longest path problem"
  [source sink graph]
  (let [s (atom {nil {:p nil :v -1}})
        order (top-order graph)
        helper (fn [s n]
                 (let [in-edges (filter #(= n (second %)) graph)
                       in-vals (map #(hash-map :p (first %)
                                               :v (+ (nth % 2) (:v (s (first %))))) in-edges)
                       v (if (not-empty in-vals) (apply max-key :v in-vals) {:p nil :v 0})]
                   (assoc s n v)))]
    (doseq [n order]
      (swap! s helper n))
    (loop [n sink path (seq [sink]) count 1000]
      (assert (pos? count))
      (if (= source n) [(apply + path) path]
        (let [p (:p (@s n))]
          (recur p (conj path p) (dec count)))))))

(defn read-graph
  [filename]
  (let [lines (split (slurp filename) #"\n")]
    (map #(let [[path, weight] (split % #":")
                [in, out] (split path #"->")]
            [(Integer/parseInt in) (Integer/parseInt out) (Integer/parseInt weight)]) lines)))

(defn blosum62
  "blosum62 scoring function"
  [a b]
  (let [potentials "ACDEFGHIKLMNPQRSTVWY"
        a-index (.indexOf potentials (str a))
        b-index (.indexOf potentials (str b))
        matrix [[ 4  0 -2 -1 -2  0 -2 -1 -1 -1 -1 -2 -1 -1 -1  1  0  0 -3 -2]
                [ 0  9 -3 -4 -2 -3 -3 -1 -3 -1 -1 -3 -3 -3 -3 -1 -1 -1 -2 -2]
                [-2 -3  6  2 -3 -1 -1 -3 -1 -4 -3  1 -1  0 -2  0 -1 -3 -4 -3]
                [-1 -4  2  5 -3 -2  0 -3  1 -3 -2  0 -1  2  0  0 -1 -2 -3 -2]
                [-2 -2 -3 -3  6 -3 -1  0 -3  0  0 -3 -4 -3 -3 -2 -2 -1  1  3]
                [ 0 -3 -1 -2 -3  6 -2 -4 -2 -4 -3  0 -2 -2 -2  0 -2 -3 -2 -3]
                [-2 -3 -1  0 -1 -2  8 -3 -1 -3 -2  1 -2  0  0 -1 -2 -3 -2  2]
                [-1 -1 -3 -3  0 -4 -3  4 -3  2  1 -3 -3 -3 -3 -2 -1  3 -3 -1]
                [-1 -3 -1  1 -3 -2 -1 -3  5 -2 -1  0 -1  1  2  0 -1 -2 -3 -2]
                [-1 -1 -4 -3  0 -4 -3  2 -2  4  2 -3 -3 -2 -2 -2 -1  1 -2 -1]
                [-1 -1 -3 -2  0 -3 -2  1 -1  2  5 -2 -2  0 -1 -1 -1  1 -1 -1]
                [-2 -3  1  0 -3  0  1 -3  0 -3 -2  6 -2  0  0  1  0 -3 -4 -2]
                [-1 -3 -1 -1 -4 -2 -2 -3 -1 -3 -2 -2  7 -1 -2 -1 -1 -2 -4 -3]
                [-1 -3  0  2 -3 -2  0 -3  1 -2  0  0 -1  5  1  0 -1 -2 -2 -1]
                [-1 -3 -2  0 -3 -2  0 -3  2 -2 -1  0 -2  1  5 -1 -1 -3 -3 -2]
                [ 1 -1  0  0 -2  0 -1 -2  0 -2 -1  1 -1  0 -1  4  1 -2 -3 -2]
                [ 0 -1 -1 -1 -2 -2 -2 -1 -1 -1 -1  0 -1 -1 -1  1  5  0 -2 -2]
                [ 0 -1 -3 -2 -1 -3 -3  3 -2  1  1 -3 -2 -2 -3 -2  0  4 -3 -1]
                [-3 -2 -4 -3  1 -2 -2 -3 -3 -2 -1 -4 -4 -2 -3 -3 -2 -3 11  2]
                [-2 -2 -3 -2  3 -3  2 -1 -2 -1 -1 -2 -3 -1 -2 -2 -2 -1  2  7]]]
    (-> matrix (nth a-index) (nth b-index))))
