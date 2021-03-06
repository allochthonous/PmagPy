FasdUAS 1.101.10   ��   ��    k             l     ��  ��    S M-----------------------------------------------------------------------------     � 	 	 � - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
  
 l     ��  ��    A ; Installation script for Lisa Tauxe's PmagPy Python Package     �   v   I n s t a l l a t i o n   s c r i p t   f o r   L i s a   T a u x e ' s   P m a g P y   P y t h o n   P a c k a g e      l     ��  ��    - ' Written by Rupert C. J. Minnett, Ph.D.     �   N   W r i t t e n   b y   R u p e r t   C .   J .   M i n n e t t ,   P h . D .      l     ��  ��      Last updated 2/10/2016     �   .   L a s t   u p d a t e d   2 / 1 0 / 2 0 1 6      l     ��  ��    S M-----------------------------------------------------------------------------     �   � - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -      l     ��������  ��  ��       !   l     �� " #��   " ( " Set the default installation path    # � $ $ D   S e t   t h e   d e f a u l t   i n s t a l l a t i o n   p a t h !  % & % l     '���� ' r      ( ) ( n      * + * 1   	 ��
�� 
psxp + l    	 ,���� , b     	 - . - l     /���� / c      0 1 0 l     2���� 2 I    �� 3��
�� .earsffdralis        afdr 3 m     ��
�� afdrcusr��  ��  ��   1 m    ��
�� 
TEXT��  ��   . m     4 4 � 5 5  P m a g P y��  ��   ) o      ���� 0 installation_path  ��  ��   &  6 7 6 l     ��������  ��  ��   7  8 9 8 l    :���� : r     ; < ; c     = > = 4    �� ?
�� 
psxf ? l    @���� @ b     A B A l    C���� C n     D E D 1    ��
�� 
psxp E l    F���� F I   �� G��
�� .earsffdralis        afdr G  f    ��  ��  ��  ��  ��   B m     H H � I I  / . .��  ��   > m    ��
�� 
ctxt < o      ���� 0 current_working_directory  ��  ��   9  J K J l     ��������  ��  ��   K  L M L l     �� N O��   N @ : Confirm that the default installation path should be used    O � P P t   C o n f i r m   t h a t   t h e   d e f a u l t   i n s t a l l a t i o n   p a t h   s h o u l d   b e   u s e d M  Q R Q l   ; S���� S r    ; T U T l   9 V���� V n    9 W X W 1   5 9��
�� 
ttxt X l   5 Y���� Y I   5�� Z [
�� .sysodlogaskr        TEXT Z l 	    \���� \ l 	    ]���� ] m      ^ ^ � _ _ � I n s t a l l   P m a g P y   t o   t h e   d e f a u l t   d i r e c t o r y   ( y o u r   h o m e   d i r e c t o r y ) ? 
 
 I f   n o t ,   p l e a s e   e d i t   t h e   i n s t a l l a t i o n   p a t h   b e l o w   a n d   p r e s s   ' O K ' .��  ��  ��  ��   [ �� ` a
�� 
appr ` l 	 ! " b���� b m   ! " c c � d d : P m a g P y   I n s t a l l a t i o n   D i r e c t o r y��  ��   a �� e f
�� 
dtxt e l 
 # $ g���� g o   # $���� 0 installation_path  ��  ��   f �� h i
�� 
btns h J   % + j j  k l k m   % & m m � n n  O K l  o�� o m   & ) p p � q q  C a n c e l��   i �� r��
�� 
dflt r m   . /���� ��  ��  ��  ��  ��   U o      ���� 0 installation_path  ��  ��   R  s t s l     ��������  ��  ��   t  u v u l     �� w x��   w 8 2 As long as the requested installation path exists    x � y y d   A s   l o n g   a s   t h e   r e q u e s t e d   i n s t a l l a t i o n   p a t h   e x i s t s v  z { z l  < � |���� | V   < � } ~ } k   I �    � � � l  I I��������  ��  ��   �  � � � l  I I�� � ���   � @ : Confirm that the default installation path should be used    � � � � t   C o n f i r m   t h a t   t h e   d e f a u l t   i n s t a l l a t i o n   p a t h   s h o u l d   b e   u s e d �  � � � r   I l � � � l  I h ����� � I  I h�� � �
�� .sysodlogaskr        TEXT � l 	 I L ����� � l 	 I L ����� � m   I L � � � � �L T h e   P m a g P y   i n s t a l l a t i o n   d i r e c t o r y   b e l o w   a l r e a d y   e x i s t s .   I s   i t   O K   t o   o v e r w r i t e   t h i s   d i r e c t o r y ? 
                 
 I f   n o t ,   p l e a s e   e d i t   t h e   i n s t a l l a t i o n   p a t h   b e l o w   a n d   p r e s s   ' O K ' .��  ��  ��  ��   � �� � �
�� 
appr � l 	 M P ����� � m   M P � � � � � H P m a g P y   I n s t a l l a t i o n   D i r e c t o r y   E x i s t s��  ��   � �� � �
�� 
dtxt � l 
 Q R ����� � o   Q R���� 0 installation_path  ��  ��   � �� � �
�� 
btns � J   S ^ � �  � � � m   S V � � � � �  O K �  � � � m   V Y � � � � �  O v e r w r i t e �  ��� � m   Y \ � � � � �  C a n c e l��   � �� ���
�� 
dflt � m   a b���� ��  ��  ��   � o      ���� 0 response   �  ��� � Z   m � � ��� � � F   m � � � � =   m x � � � n   m t � � � 1   p t��
�� 
bhit � o   m p���� 0 response   � m   t w � � � � �  O v e r w r i t e � =   { � � � � n   { � � � � 1   ~ ���
�� 
ttxt � o   { ~���� 0 response   � o   � ����� 0 installation_path   � I  � ��� ���
�� .sysoexecTEXT���     TEXT � b   � � � � � m   � � � � � � �  r m   - R f   � n   � � � � � 1   � ���
�� 
strq � o   � ����� 0 installation_path  ��  ��   � r   � � � � � n   � � � � � 1   � ���
�� 
ttxt � o   � ����� 0 response   � o      ���� 0 installation_path  ��   ~ I   @ H�� ����� 0 
fileexists 
fileExists �  ��� � n   A D � � � 1   B D��
�� 
psxp � o   A B���� 0 installation_path  ��  ��  ��  ��   {  � � � l     ��������  ��  ��   �  � � � l     �� � ���   � Z T The installation directory either didn't exist or has been deleted, so create a one    � � � � �   T h e   i n s t a l l a t i o n   d i r e c t o r y   e i t h e r   d i d n ' t   e x i s t   o r   h a s   b e e n   d e l e t e d ,   s o   c r e a t e   a   o n e �  � � � l  � � ����� � I  � ��� ���
�� .sysoexecTEXT���     TEXT � b   � � � � � m   � � � � � � �  m k d i r   - p   � l  � � ���� � n   � � � � � 1   � ��~
�~ 
strq � n   � � � � � 1   � ��}
�} 
psxp � o   � ��|�| 0 installation_path  ��  �  ��  ��  ��   �  � � � l     �{�z�y�{  �z  �y   �  � � � l     �x � ��x   � L F Copy the current directory's contents into the installation directory    � � � � �   C o p y   t h e   c u r r e n t   d i r e c t o r y ' s   c o n t e n t s   i n t o   t h e   i n s t a l l a t i o n   d i r e c t o r y �  � � � l  � � ��w�v � I  � ��u ��t
�u .sysoexecTEXT���     TEXT � b   � � � � � b   � � � � � b   � � � � � m   � � � � � � �  c p   - R   � l  � � ��s�r � n   � � � � � 1   � ��q
�q 
strq � n   � � � � � 1   � ��p
�p 
psxp � o   � ��o�o 0 current_working_directory  �s  �r   � m   � � � � � � �    � l  � � ��n�m � n   � � � � � 1   � ��l
�l 
strq � n   � � � � � 1   � ��k
�k 
psxp � o   � ��j�j 0 installation_path  �n  �m  �t  �w  �v   �  � � � l     �i�h�g�i  �h  �g   �  � � � l     �f � �f   � W Q Permanently append the installation path to the user's PATH environment variable     � �   P e r m a n e n t l y   a p p e n d   t h e   i n s t a l l a t i o n   p a t h   t o   t h e   u s e r ' s   P A T H   e n v i r o n m e n t   v a r i a b l e �  l  � ��e�d I  � ��c
�c .sysoexecTEXT���     TEXT b   � � b   � �	
	 b   � � b   � � m   � � � 
 e c h o   l  � ��b�a n   � � 1   � ��`
�` 
strq n   � � 1   � ��_
�_ 
psxp o   � ��^�^ 0 installation_path  �b  �a   m   � � �  :
 l  � ��]�\ b   � � n   � � 1   � ��[
�[ 
strq n   � � 1   � ��Z
�Z 
psxp o   � ��Y�Y 0 installation_path   m   � � �    / p r o g r a m s�]  �\   m   � �!! �"" ,   >   / e t c / p a t h s . d / P m a g P y �X#�W
�X 
badm# m   � ��V
�V boovtrue�W  �e  �d   $%$ l     �U�T�S�U  �T  �S  % &'& l     �R()�R  ( ] W Permanently append the installation path to the user's PYTHONPATH environment variable   ) �** �   P e r m a n e n t l y   a p p e n d   t h e   i n s t a l l a t i o n   p a t h   t o   t h e   u s e r ' s   P Y T H O N P A T H   e n v i r o n m e n t   v a r i a b l e' +,+ l     �Q-.�Q  - � �do shell script "launchctl unsetenv $PYTHONPATH:" & (quoted form of POSIX path of installation_path) with administrator privileges   . �// d o   s h e l l   s c r i p t   " l a u n c h c t l   u n s e t e n v   $ P Y T H O N P A T H : "   &   ( q u o t e d   f o r m   o f   P O S I X   p a t h   o f   i n s t a l l a t i o n _ p a t h )   w i t h   a d m i n i s t r a t o r   p r i v i l e g e s, 010 l     �P�O�N�P  �O  �N  1 232 l     �M45�M  4 ) # Ask if a PmagPy test should be run   5 �66 F   A s k   i f   a   P m a g P y   t e s t   s h o u l d   b e   r u n3 787 l  �9�L�K9 r   �:;: l  �<�J�I< I  ��H=>
�H .sysodlogaskr        TEXT= l 	 �?�G�F? l 	 �@�E�D@ m   �AA �BB P m a g P y   i n s t a l l e d   s u c c e s s f u l l y ! 
 
 W o u l d   y o u   l i k e   t o   r u n   a   P m a g P y   t e s t   ( e q u i v a l e n t   t o   e x e c u t i n g   ' e q a r e a . p y   - h '   o n   t h e   c o m m a n d   l i n e ) ?�E  �D  �G  �F  > �CCD
�C 
apprC l 	E�B�AE m  FF �GG : P m a g P y   I n s t a l l e d   S u c c e s s f u l l y�B  �A  D �@HI
�@ 
btnsH J  JJ KLK m  	MM �NN  Y e sL O�?O m  	PP �QQ  N o�?  I �>R�=
�> 
dfltR m  �<�< �=  �J  �I  ; o      �;�; 0 response  �L  �K  8 STS l CU�:�9U Z  CVW�8�7V =  (XYX n  $Z[Z 1   $�6
�6 
bhit[ o   �5�5 0 response  Y m  $'\\ �]]  Y e sW O  +?^_^ k  1>`` aba I 16�4�3�2
�4 .miscactvnull��� ��� null�3  �2  b c�1c I 7>�0d�/
�0 .coredoscnull��� ��� ctxtd m  7:ee �ff  e q a r e a . p y   - h�/  �1  _ m  +.gg�                                                                                      @ alis    l  Macintosh HD               ���DH+  �Terminal.app                                                   Ck����        ����  	                	Utilities     ���      ���"    �N  2Macintosh HD:Applications: Utilities: Terminal.app    T e r m i n a l . a p p    M a c i n t o s h   H D  #Applications/Utilities/Terminal.app   / ��  �8  �7  �:  �9  T hih l     �.�-�,�.  �-  �,  i jkj l DFl�+�*l L  DF�)�)  �+  �*  k mnm l     �(�'�&�(  �'  �&  n o�%o i     pqp I      �$r�#�$ 0 
fileexists 
fileExistsr s�"s o      �!�! 0 	posixfile 	posixFile�"  �#  q L     tt H     uu c     vwv l    x� �x c     yzy l    {��{ I    �|�
� .sysoexecTEXT���     TEXT| b     }~} b     � m     �� ���  t e s t   - e  � n    ��� 1    �
� 
strq� o    �� 0 	posixfile 	posixFile~ m    �� ���    ; e c h o   $ ?�  �  �  z m    �
� 
long�   �  w m    �
� 
bool�%       
�����������  � ��������� 0 
fileexists 
fileExists
� .aevtoappnull  �   � ****� 0 installation_path  � 0 current_working_directory  � 0 response  �  �  �  � �
q�	�����
 0 
fileexists 
fileExists�	 ��� �  �� 0 	posixfile 	posixFile�  � �� 0 	posixfile 	posixFile� ������ 
� 
strq
� .sysoexecTEXT���     TEXT
� 
long
�  
bool� ��,%�%j �&�&� �����������
�� .aevtoappnull  �   � ****� k    F��  %��  8��  Q��  z��  ���  ��� �� 7�� S�� j����  ��  ��  �  � 4������ 4������ H���� ^�� c���� m p���������� � � � � ����� ��� ����� � � �!��AFMP��\g��e��
�� afdrcusr
�� .earsffdralis        afdr
�� 
TEXT
�� 
psxp�� 0 installation_path  
�� 
psxf
�� 
ctxt�� 0 current_working_directory  
�� 
appr
�� 
dtxt
�� 
btns
�� 
dflt�� 
�� .sysodlogaskr        TEXT
�� 
ttxt�� 0 
fileexists 
fileExists�� 0 response  
�� 
bhit
�� 
bool
�� 
strq
�� .sysoexecTEXT���     TEXT
�� 
badm�� 
�� .miscactvnull��� ��� null
�� .coredoscnull��� ��� ctxt��G�j �&�%�,E�O*�)j �,�%/�&E�O�������a lva ka  a ,E�O lh*��,k+ a �a ���a a a mva ka  E` O_ a ,a  	 _ a ,� a & a �a  ,%j !Y _ a ,E�[OY��Oa "��,a  ,%j !Oa #��,a  ,%a $%��,a  ,%j !Oa %��,a  ,%a &%��,a  ,a '%%a (%a )el !Oa *�a +�a ,a -lva ka . E` O_ a ,a /  a 0 *j 1Oa 2j 3UY hOh� ��� ( / U s e r s / c r o w a n / P m a g P y� ��� L M a c i n t o s h   H D : U s e r s : c r o w a n : C o d e : P m a g P y :� ����
�� 
bhit� ���  O v e r w r i t e� �����
�� 
ttxt� ��� ( / U s e r s / c r o w a n / P m a g P y��  �  �  �   ascr  ��ޭ