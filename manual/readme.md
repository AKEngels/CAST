# CAST Manual

Verified compilation of the manual on Windows 10 using MikTex (http://miktex.org/download) with TeXStudio-Environment (http://www.texstudio.org/).

Please note than different versions of the manual may be rendered through the use of conditionals, namely:

- Verbose (More content)
- Development (Notes on the code base, compilation and more information suitable for people with access to the source code)
- Dev-Mode (Includes notes on unfinished paragraphes etc)

The extra sections enclosed by these conditionals may be activated or deactived through:

\devmodetrue
\verbosetrue
\developmenttrue

or, respectivly,

\devmodefalse
\verbosefalse
\developmentfalse

If you write notes on other developers into the Dev-Mode, it is a good idea to highlight these comments thorugh the use of colors, for example: \colorbox{red}{here_be_content}