;; This script lints all .f90 files in the directories [./src/, ./app/, and ./test/].
;; It applies the default linting rules of the emacs f90-mode.
;;
;; It is recommended that you run the linter before commiting your changes.
;;
;; Run the linter by saying the following in your shell: emacs --batch -l linter.el
;;
;; Note that if you have a currently open buffer visiting the file that has been linted, you will not be able to see the changes unless you "refresh" the view. To refresh, say C-x C-v on the file.

(defun lint-all-files (directory)
  (let ((file-list (directory-files directory t "\\.f90$")))
    (dolist (file file-list)
      (with-temp-buffer
        (insert-file-contents file)
        (f90-mode)
        (indent-region (point-min) (point-max))  ;; "TAB" entire buffer
        (write-file file)))))

(mapcar #'lint-all-files '("./src/" "./app/" "./test/" "./V3gpu/"))
