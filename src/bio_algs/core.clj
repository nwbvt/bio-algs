(ns bio-algs.core
  (:require [clojure.java.io :refer :all])
  (:require [clojure.string :refer [split trim join]]))

(defn write-result
  "Writes the solution in the expected format"
  [values]
  (with-open [outfile (writer "result.txt")]
    (binding [*out* outfile]
      (print (join "\n" values)))))

(defn read-strands
  "Reads a list of \n deliminated dna strands"
  [input]
  (map trim (split input #"\n")))
