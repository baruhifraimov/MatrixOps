# מטלה מספר 2 - מטריצות - העמסת אופרטורים - הגשה ביחידים

עבר הרבה זמן מאז שסיימתם את הקורס אלגברה לינארית :-) 
כדי לעורר בכם נוסטלגיה, במטלה הזאת אתם תממשו מחלקה המייצגת **מטריצות ריבועיות**.

כתבו מחלקה בשם ``SquareMat`` המייצגת מטריצה ריבועית. עליכם לממש במחלקה את הפונקציות הבאות:

- אופרטורי חיבור וחיסור (`+`,`-`) - מבצע חיסור איברים בהתאם למקומם בין שתי מטריצות בעלות אותו גודל.
- אופרטור מינוס אונארי (`-`) - הופך את סימן כל האיברים במטריצה (חיובי לשלילי ולהפך).
- אופרטורי השמה משולבים (`+=`, `-=`, `*=`, `/=`, `%=`)  - מבצע פעולות חיבור, חיסור, כפל, חלוקה או מודולו על המטריצה ומעדכן אותה ישירות.
- אופרטור כפל (`*`) - מבצע כפל בין 2 מטריצות.
- כפל איברי (`%`) - מכפיל כל איבר במטריצה אחת באיבר המתאים במטריצה השנייה.
- חילוק בסקלר (`/`) - מחלק כל איבר במטריצה במספר סקלרי.
- מודולו עם סקלר (`% int`) - מבצע פעולה של מודולו על כל איבר במטריצה עם מספר שלם נתון.
- אופרטור חזקה (`^`) - מעלה את המטריצה בחזקה על ידי כפל חוזר של המטריצה בעצמה.
- אופרטורי הגדלה והקטנה ב-1 (`++`, `--`) - מגדילים או מקטינים את כל איברי המטריצה באחד.
- אופרטור Transpose (`~`) - מחליף את השורות בעמודות של המטריצה.
- גישה לאיברים באמצעות שני אינדקסים (`()`) - מאפשר גישה ישירה לאיבר ספציפי במטריצה על ידי שימוש באינדקסים `(i, j)`.
- אופרטורי שוויון ו-אי-שוויון (`==`,`!=`) - ככלל - מטריצות הן שוות אם סכום האיברים שלהם זהה.
- אופרטורי השוואה (`<`, `>`, `<=`, `>=`) - משווה בין שתי מטריצות לפי כלל מסוים, למשל סכום האיברים שלהן. ככלל - מטריצה גדולה מהשנייה אם הסכום של האיברים שלה גדול מסכום האיברים של המטריצה השנייה.
- אופרטור דטרמיננטה - (`!`) - מחשב את הדטרמיננטה של המטריצה.
- אופרטור פלט - (`>>`) - מדפיס את המטריצה בצורה הגיונית.

יש לבדוק את תקינות הקלט ולזרוק שגיאות במידת הצורך. בנוסף עליכם לכתוב טסטים לכל הפונקציות שכתבתם.

מצורף קובץ בשם ``Demo.cpp`` המדגים את הקוד.

יש לכתוב קובץ main שבו אתם מדגימים את הקוד שכתבתם. בנוסף, עליכם ליצור קובץ ``Makefile`` כאשר הפקודה ``make Main`` 
תריץ את הקובץ לדוגמה שכתבתם והפקודה ``make test`` תריץ את הטסטים שלכם.

כמו כן, עליכם לכתוב בתחילת **כל** קובץ את המייל שלכם. עליכם להגיש קובץ ``README`` המסביר את ההיררכיה של המחלקות ובאילו ספריות השתמשתם. כמו כן, עליכם לכתוב בתחילת כל קובץ את המייל שלכם.

כדי להגיש את המטלה במודל עליכם להגיש קובץ טקסט המורכב מ-3 שורות: בשורה הראשונה מספר תעודת הזהות שלכם, בשורה השנייה הקישור להגשה שלכם ובשורה השלישית הקומיט האחרון של הפרויקט.

מותר להתייעץ עם סטודנטים אחרים או לחפש מידע באינטרנט; אסור להעתיק קטעי-קוד שלמים מסטודנטים אחרים, מהאינטרנט וממודלי בינה מלאכותית (למשל, chatGPT).
יש לדווח על כל עזרה שקיבלתם, מסטודנטים אחרים או מהאינטרנט, בהתאם ל[תקנון היושר של המחלקה](https://www.ariel.ac.il/wp/cs/wp-content/uploads/sites/88/2020/08/Guidelines-for-Academic-Integrity.pdf).

בהצלחה!
