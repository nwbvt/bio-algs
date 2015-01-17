(ns bio-algs.rearrangements
  (:require [clojure.string :refer [join]]))

(defn reversal
  "Do a reversal of a permutation"
  [perm start end]
  (assert (not (neg? start)))
  (assert (not (neg? end)))
  (let [[left postfix] (split-at end perm)
        [prefix inner] (split-at start left)
        reversed (map - (reverse inner))]
    (concat prefix reversed postfix)))

(defn greedy-sorting
  "Given a permutation, returns the sequence of permutations corresponding to the greedy sorting algorithm"
  [initial]
  (loop [perm initial perms [] i 1]
    (if (> i (count initial)) perms
      (if (= i (nth perm (dec i)))
        (recur perm perms (inc i))
        (let [index (max (.indexOf perm i) (.indexOf perm (- i))) ; the one that doesn't exist will be -1
              _ (assert (not (neg? index)))
              start (dec i)
              end (inc index)
              next-perm (reversal perm start end)]
          (recur next-perm (conj perms next-perm) i))))))

(defn format-perm
  "Formats a permutation"
  [perm]
  (str "(" (join " " (map #(format "%+d" %)perm)) ")"))

(defn count-breakpoints
  "counts the number of breakpoints in the permutation"
  [perm]
  (let [n (count perm)
        with-edges (concat [0] perm [(inc n)])]
    (count (filter #(not= 1 %) (map #(- (nth with-edges (inc %)) (nth with-edges %)) (range (inc n)))))))

(defn find-connected-vertex
  "finds the connected vertex for a given vertex in a given genome
   the - vertex is the vertex incoming to a given block, + outgoing"
  [vertex genome]
  (some (fn [chromosome]
          (let [indices (map #(.indexOf chromosome %) [vertex (- vertex)])
                [next-index, mult] (cond
                                     (<= 0 (first indices)) [(inc (first indices)), -1]
                                     (<= 0 (second indices)) [(dec (second indices)), 1])]
            (if next-index (* mult (nth chromosome (mod next-index (count chromosome)))))))
        genome))

(defn- rem-cycles
  "recusion function for get cycles"
  [vert left genomes cur]
  (if (empty? left) [cur]
    (let [next-vert (find-connected-vertex vert (first genomes))]
      (if (left next-vert)
        (rem-cycles next-vert (disj left next-vert) (rest genomes) (conj cur next-vert))
        (cons cur (lazy-seq (rem-cycles (first left) (set (rest left)) (rest genomes) [(first left)])))))))

(defn get-cycles
  "gets the cycles in the breakpoint graph between two genomes"
  [genome1 genome2]
  (let [blocks (apply concat genome1)
        all-verts (set (concat blocks (map - blocks)))]
    (rem-cycles (first all-verts) (set (rest all-verts)) (cycle [genome1 genome2]) [(first all-verts)])))

(defn count-cycles
  "counts the number of cycles in the breakpoint graph between two genomes"
  [genome1 genome2]
  (count (get-cycles genome1 genome2)))


(defn distance
  "Calculate the 2 break distance between two genomes"
  [genome1 genome2]
  (- (apply + (map count genome1))
     (count-cycles genome1 genome2)))

(defn non-trivial?
  "True if the given cycle is not trivial"
  [cyc]
  (< 2 (count cyc)))

(defn make-break
  "Makes a 2 break change to a genome"
  [genome break1 break2]
  )

(defn two-break-sort
  "Sorts using two break sorting"
  [start end]
  (assert (= (->> start (apply concat) (map #(Math/abs %)) sort)
             (->> end (apply concat) (map #(Math/abs %)) sort))
          "Genomes do not contain same blocks, cannot be sorted")
  (loop [genome start order [start]]
    (let [cyc (first (filter non-trivial? (get-cycles genome end)))]
      (if (nil? cyc) order
        
        )
      )))
