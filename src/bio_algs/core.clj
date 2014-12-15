(ns bio-algs.core
  (:require [clojure.java.io :refer :all])
  (:require [clojure.string :refer [split trim join]]))

(defn write-result
  "Writes the solution in the expected format"
  [values]
  (with-open [outfile (writer "result.txt")]
    (binding [*out* outfile]
      (print (join "\n" values)))))

(defn write-1-result
  "Writes a single result in the expected format"
  [value]
  (write-result [value]))

(defn read-strands
  "Reads a list of \n deliminated dna strands"
  [input]
  (map trim (split input #"\n")))

(defn read-file
  "Reads strings from a file"
  [filename]
  (read-strands (apply str (slurp filename))))
