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

(defn get-start
  "gets the starting position for the backtrack"
  [s max-i max-j locality]
  (case locality
    :local (apply max-key (fn [[a b]] (:v (s [a b]))) (for [a (range max-i) b (range max-j)] [a b])) 
    :global [max-i max-j] 
    :fit [(apply max-key #(:v (s [% max-j])) (range max-i)) max-j]))

(defn find-backtrack
  "backtracks through the graph to find the longest common subseq"
  [s seq1 seq2 locality]
  (let [[start-i start-j] (get-start s (count seq1) (count seq2) locality)]
   (loop [i start-i j start-j align1 '() align2 '()]
    (if (and (zero? i) (zero? j)) [(:v (s [start-i start-j])) (apply str align1) (apply str align2)]
      (let [v (cond (zero? i) (if (= :local locality) {:d :start} {:d :right})
                    (zero? j) (if (not= :global locality) {:d :start} {:d :down})
                    :else (s [i, j]))]
        (case (:d v)
          :start  (recur 0 0 align1 align2)
          :down   (recur (dec i) j (conj align1 (nth seq1 (dec i))) (conj align2 "-"))
          :right  (recur i (dec j) (conj align1 "-") (conj align2 (nth seq2 (dec j))))
          :across (recur (dec i) (dec j) (conj align1 (nth seq1 (dec i))) (conj align2 (nth seq2 (dec j))))))))))

(defn simple-matching
  "simple matching"
  [penalty i j]
  (if (= i j) 1 (* -1 penalty)))

(defn longest-common-subseq
  "find the longest common subsequence of two dna strands"
  ([seq1 seq2] (longest-common-subseq seq1 seq2 0))
  ([seq1 seq2 indel-pen] (longest-common-subseq seq1 seq2 indel-pen (partial simple-matching 0)))
  ([seq1 seq2 indel-pen score-fn] (longest-common-subseq seq1 seq2 indel-pen score-fn :global))
  ([seq1 seq2 indel-pen score-fn locality]
    (let [s (atom {[0 0] {:v 0}})
          helper (fn [s i j]
                   (assoc s [i,j]
                     (let [start {:v (if (or (= :local locality) (and (= :fit locality) (= 0 j))) 0 Double/NEGATIVE_INFINITY) :d :start}
                           d-val {:v (if (zero? i) Double/NEGATIVE_INFINITY (- (:v (s [(dec i) j])) indel-pen)) :d :down}
                           r-val {:v (if (zero? j) Double/NEGATIVE_INFINITY (- (:v (s [i (dec j)])) indel-pen)) :d :right}
                           across-val {:v (if (or (zero? i) (zero? j)) Double/NEGATIVE_INFINITY
                                            (+ (score-fn (nth seq1 (dec i)) (nth seq2 (dec j))) 
                                               (:v (s [(dec i) (dec j)]))))
                                       :d :across}]
                       (max-key :v d-val r-val across-val start))))]
      (doseq [i (range 0 (inc (count seq1))) j (range 0 (inc (count seq2))) :when (or (pos? i) (pos? j))]
        (swap! s helper i j))
      (find-backtrack @s seq1 seq2 locality))))

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

(defn edit-distance
  "Finds the edit distance between 2  strings"
  [seq1 seq2]
  (* -1 (first (longest-common-subseq seq1 seq2 1 (fn [a b] (if (= a b) 0 -1))))))

(defn score-mat
  "Creates a scoring function from a matrix"
  [matrix]
  (fn [a b]
    (let [potentials "ACDEFGHIKLMNPQRSTVWY"
          a-index (.indexOf potentials (str a))
          b-index (.indexOf potentials (str b))]
      (-> matrix (nth a-index) (nth b-index)))))

(def blosum62
  (score-mat [[ 4  0 -2 -1 -2  0 -2 -1 -1 -1 -1 -2 -1 -1 -1  1  0  0 -3 -2]
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
              [-2 -2 -3 -2  3 -3  2 -1 -2 -1 -1 -2 -3 -1 -2 -2 -2 -1  2  7]]))

(def pam250
  (score-mat [[ 2 -2  0  0 -3  1 -1 -1 -1 -2 -1  0  1  0 -2  1  1  0 -6 -3]
              [-2 12 -5 -5 -4 -3 -3 -2 -5 -6 -5 -4 -3 -5 -4  0 -2 -2 -8  0]
              [ 0 -5  4  3 -6  1  1 -2  0 -4 -3  2 -1  2 -1  0  0 -2 -7 -4]
              [ 0 -5  3  4 -5  0  1 -2  0 -3 -2  1 -1  2 -1  0  0 -2 -7 -4]
              [-3 -4 -6 -5  9 -5 -2  1 -5  2  0 -3 -5 -5 -4 -3 -3 -1  0  7]
              [ 1 -3  1  0 -5  5 -2 -3 -2 -4 -3  0  0 -1 -3  1  0 -1 -7 -5]
              [-1 -3  1  1 -2 -2  6 -2  0 -2 -2  2  0  3  2 -1 -1 -2 -3  0]
              [-1 -2 -2 -2  1 -3 -2  5 -2  2  2 -2 -2 -2 -2 -1  0  4 -5 -1]
              [-1 -5  0  0 -5 -2  0 -2  5 -3  0  1 -1  1  3  0  0 -2 -3 -4]
              [-2 -6 -4 -3  2 -4 -2  2 -3  6  4 -3 -3 -2 -3 -3 -2  2 -2 -1]
              [-1 -5 -3 -2  0 -3 -2  2  0  4  6 -2 -2 -1  0 -2 -1  2 -4 -2]
              [ 0 -4  2  1 -3  0  2 -2  1 -3 -2  2  0  1  0  1  0 -2 -4 -2]
              [ 1 -3 -1 -1 -5  0  0 -2 -1 -3 -2  0  6  0  0  1  0 -1 -6 -5]
              [ 0 -5  2  2 -5 -1  3 -2  1 -2 -1  1  0  4  1 -1 -1 -2 -5 -4]
              [-2 -4 -1 -1 -4 -3  2 -2  3 -3  0  0  0  1  6  0 -1 -2  2 -4]
              [ 1  0  0  0 -3  1 -1 -1  0 -3 -2  1  1 -1  0  2  1 -1 -2 -3]
              [ 1 -2  0  0 -3  0 -1  0  0 -2 -1  0  0 -1 -1  1  3  0 -5 -3]
              [ 0 -2 -2 -2 -1 -1 -2  4 -2  2  2 -2 -1 -2 -2 -1  0  4 -6 -2]
              [-6 -8 -7 -7  0 -7 -3 -5 -3 -2 -4 -4 -6 -5  2 -2 -5 -6 17  0]
              [-3  0 -4 -4  7 -5  0 -1 -4 -1 -2 -2 -5 -4 -4 -3 -3 -2  0 10]]))
