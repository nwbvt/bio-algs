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

(defn find-next-block
  "finds the next block for a given block in a given genome"
  [block genome]
  (some (fn [chromosome]
          (let [indices (map #(.indexOf chromosome %) [block (- block)])
                [next-index, mult] (cond
                                     (<= 0 (first indices)) [(inc (first indices)), 1]
                                     (<= 0 (second indices)) [(dec (second indices)), -1])]
            (if next-index (* mult (nth chromosome (mod next-index (count chromosome)))))))
        genome))

(defn count-cycles
  "counts the number of cycles in the breakpoint graph between two genomes"
  [genome1 genome2]
  (let [blocks (apply concat genome1)
        all-verts (set (concat blocks (map - blocks)))]
    (loop [cycles 0 vert (first all-verts) left (set (rest all-verts)) genomes (cycle [genome1 (map reverse genome2)])]
      (if (empty? left) cycles
        (let [next-vert (find-next-block vert (first genomes))]
          (if (left next-vert)
            (recur cycles next-vert (disj left next-vert) (rest genomes))
            (recur (inc cycles) (first left) (set (rest left)) (rest genomes))))))))

(defn two-break
  "Calculate the 2 break distance between two genomes"
  [genome1 genome2]
  (- (apply + (map count genome1))
     (count-cycles genome1 genome2)))
