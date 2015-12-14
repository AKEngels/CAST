# cast (Conformational Analysis and Search Tool)
This is the official repository for the
CAST
Conformational Analysis and Search Tool.

CAST currently has no external dependencies and reads its runtime options from the INPUTFILE. CAST's features are split into tasks that are outlined in the manual (hopefully).

Notes for (new) developers:

On coding CAST:
- Try to make your code understanible, use variable names that reflect the usage and purpose and so on.

On using the git repositry:
- Make commits atomic (as small as possible). This helps tracking bugs.
- Write useful commit messages (see commit history for examples)
- If your desired commit breaks backwards compatability or simply does not compile: COMMIT TO A BRANCH!
- The master branch should be the current "nightly" build and should always compile.
- Use issues and also resolve some issues if you have time. :)
- Don't hesitate to comment on commits or issues if you have questions or suggestions
