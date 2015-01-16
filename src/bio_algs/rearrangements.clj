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
