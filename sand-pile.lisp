;;; sand-pile.lisp --- 2D CA simulating a sand pile

;; Copyright 2012 Eric Schulte

;; This program is free software: you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or
;; (at your option) any later version.

;; This program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.

;; You should have received a copy of the GNU General Public License
;; along with this program.  If not, see <http://www.gnu.org/licenses/>.

;;; Commentary:

;; This program implements the 2D CA used to investigate
;; self-organized criticality through the simulation of a simplified
;; sand-pile model in the following publication.

;; @article{bak1988self,
;;   author={Bak, P. and Tang, C. and Wiesenfeld, K. and others},
;;   title={Self-organized criticality},
;;   journal={Physical review A},
;;   year={1988},
;;   volume={38},
;;   number={1},
;;   pages={364--374}
;; }

;;; Code:
(require 'alexandria) (require 'metabang-bind)
(defpackage #:sand-pile (:use :common-lisp :alexandria :metabang-bind))
(in-package :sand-pile)

(defvar max-init       10)
(defvar critical-slope 2)
(defvar *size*         20)
(defvar *dim*          2)
(defvar *grid*
  (make-array (make-list *dim* :initial-element *size*) :element-type :number))

(defun z (coord) (apply #'aref *grid* coord))
(defun set-z (coord new) (setf (apply #'aref *grid* coord) new))
(defsetf z set-z)

(defun init ()
  "Randomly initialize the `*grid*'."
  (mapc (lambda (coord) (setf (z coord) (random max-init)))
        (cross (array-dimensions *grid*))))

(defun update-rule (coord)
  "Update a coordinate in the `*grid*'."
  (mapc (lambda (neighbor)
          (when (> (- (z coord) (z neighbor)) critical-slope)
            (decf (z coord))
            (incf (z neighbor))))
        (sort (neighbors coord) #'< :key (lambda (it) (apply #'aref *grid* it)))))

(defun run (steps)
  (init)
  (let ((history (list (copy-array *grid*))))
    (dotimes (n steps history)
      (mapc #'update-rule (shuffle (cross (array-dimensions *grid*))))
      (incf (z (random-elt (cross (array-dimensions *grid*)))))
      (push (copy-array *grid*) history))))


;;; Utility
(defun range (n)
  (loop for i upto (1- n) collect i))

(defun cross (dimensions)
  (if (cdr dimensions)
      (let ((this (range (car dimensions))))
        (mapcan (lambda (rst) (mapcar (lambda (el) (cons el rst)) this))
                (cross (cdr dimensions))))
      (mapcar #'list (range (car dimensions)))))

(defun neighbors (dim)
  (remove-if (lambda (dim) (some (lambda (n) (or (< n 0) (<= *size* n))) dim))
             (mapcar (lambda (off) (map 'list #'+ off dim))
                     (mapcar (lambda (off) (mapcar #'1- off))
                             (remove (make-list *dim* :initial-element 1)
                               (cross (make-list *dim* :initial-element 3))
                               :test #'tree-equal)))))

;;; sand-pile.lisp ends here
