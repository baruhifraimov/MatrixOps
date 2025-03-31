<div dir="rtl">

# מטלה מספר 2 - מימוש אופרטורים עבור מטריצות ריבועיות

### יושרה אקדמית

במהלך העבודה על המטלות, מותר להתייעץ עם סטודנטים אחרים ולחפש מידע באינטרנט. עם זאת, חל איסור להעתיק קטעי קוד שלמים ממקורות חיצוניים, כולל סטודנטים אחרים, אתרי אינטרנט ומודלי בינה מלאכותית (כגון ChatGPT).

יש לדווח על כל עזרה שקיבלתם, בין אם מדובר בהתייעצות עם סטודנטים אחרים או במידע שנמצא באינטרנט, בהתאם ל[תקנון היושר של המחלקה](https://www.ariel.ac.il/wp/cs/wp-content/uploads/sites/88/2020/08/Guidelines-for-Academic-Integrity.pdf).
**במקרה של שימוש בכלי בינה מלאכותית (AI), יש לצרף את הפרומפטים שהוזנו ואת התשובות שהתקבלו**.

-----
* **מטרת המטלה:** הבנת החומר הנלמד בהרצאות השלישית והרביעית, כגון: העמסת אופרטורים ופונקציות חברות.
* **שימו לב!** ההגשה של המטלה ביחידים.

---

## הוראות הגשה ב Moodle:

במערכת Moodle יש להגיש **קובץ טקסט למשל (`submission.txt`)** המכיל 3 שורות בפורמט הבא:

1. **תעודת זהות** – מספר תעודת הזהות של הסטודנט.
2. **קישור להגשה** – קישור למאגר ה-GitHub שבו נמצא הפרויקט.
3. **פרטי ה-commit האחרון** – המחרוזת המזהה של ה-commit האחרון (`commit hash`) 

 - דוגמה לקובץ הגשה תקין:
```
123456789
https://github.com/example-user/assignment
e3f1c1a 
```

---

עבר הרבה זמן מאז שסיימתם את הקורס אלגברה לינארית :-) 
כדי לעורר בכם נוסטלגיה, במטלה הזאת אתם תממשו מחלקה המייצגת **מטריצות ריבועיות** של מספרים ממשיים.

---

## דרישות המטלה:

### מימוש המחלקות הבאות:

הוסיפו את המחלקות במרחב שמות (namespace) חדש לבחירתכם.

#### מחלקה בשם `SquareMat`:
המטריצה היא מטריצה ריבועית המכילה מספרים ממשיים. עבור קלט לא תקין יש לזרוק חריגה.
המחלקה תכיל את המטריצה וכן את הפונקציות הבאות:

- אופרטורי חיבור וחיסור (`+`,`-`) - מבצעים חיבור/חיסור איברים בהתאם למקומם בין שתי מטריצות בעלות אותו גודל.
- אופרטור מינוס אונארי (`-`) - הופך את סימן כל האיברים במטריצה (חיובי לשלילי ולהפך).
- אופרטורי השמה משולבים (`=+`, `=-`, `=*`, `=/`, `=%`)  - מבצעים פעולות של חיבור, חיסור, כפל, חלוקה או מודולו על המטריצה ומעדכנים אותה ישירות.
- אופרטור כפל (`*`) - מבצע [כפל בין 2 מטריצות](https://he.wikipedia.org/wiki/%D7%9B%D7%A4%D7%9C_%D7%9E%D7%98%D7%A8%D7%99%D7%A6%D7%95%D7%AA).
- אופרטור כפל בסקלר (`scalar*matrix` וגם `matrix*scalar`)
- כפל איברי (`%`) - מכפיל כל איבר במטריצה אחת באיבר המתאים במטריצה השנייה.
- חילוק בסקלר (`/`) - מחלק כל איבר במטריצה במספר סקלרי.
- מודולו עם סקלר (`%`) - מבצע פעולה של מודולו על כל איבר במטריצה עם מספר שלם נתון.
- אופרטור חזקה (`^`) - מעלה את המטריצה בחזקה על ידי כפל חוזר של המטריצה בעצמה.
- אופרטורי הגדלה והקטנה ב-1 (`++`, `--`) - מגדילים או מקטינים את כל איברי המטריצה באחד.
- אופרטור Transpose (`~`) - מחליף את השורות בעמודות של המטריצה.
- גישה לאיברים באמצעות שני אינדקסים (`()`) - מאפשר גישה ישירה לאיבר ספציפי במטריצה על ידי שימוש באינדקסים `(i, j)`.
- אופרטורי שוויון ו-אי-שוויון (`==`,`=!`) - ככלל - מטריצות הן שוות אם סכום האיברים שלהם זהה.
- אופרטורי השוואה (`<`, `>`, `=>`, `=<`) - משווה בין שתי מטריצות לפי כלל מסוים, למשל סכום האיברים שלהן. ככלל - מטריצה גדולה מהשנייה אם הסכום של האיברים שלה גדול מסכום האיברים של המטריצה השנייה.
- אופרטור דטרמיננטה - (`!`) - מחשב את הדטרמיננטה של המטריצה.
- אופרטור פלט - (`>>`) - מדפיס את המטריצה בצורה הגיונית.

* יש לבדוק את תקינות הקלט ולזרוק שגיאות במידת הצורך. בנוסף עליכם לכתוב טסטים לכל הפונקציות שכתבתם.


---


#### דרישות נוספות:

- חשוב לוודא שה-repository ציבורי.
- כתבו בתחילת **כל** קובץ את כתובת המייל שלכם.
- כתבו קוד נקי, מסודר, מחולק לקבצים, מודולרי, מתועד בצורה מספקת וכמובן בדיקות יחידה עבור כל הפונקציות.
- בדקו את תקינות הקלט ולזרוק חריגות מתאימות במידת הצורך.
- לשימושכם הקישור הבא [doctest](https://github.com/doctest/doctest) בו תוכלו לראות דוגמאות נוספות לשימוש בסיפריה זו.
- יש לבדוק שאין זליגת זיכרון באמצעות `valgrind`.
- יש לצרף גם קובץ `README` עם הסבר על פרויקט, על חלוקה למחלקות וקבצים וכל מידע אחר רלוונטי.


#### קובץ `Makefile`:
הוסיפו לפרויקט קובץ `Makefile` הכולל את הפקודות הבאות:
- הפקודה `make Main` – להרצת קובץ ההדגמה.
- הפקודה `make test` – להרצת בדיקות היחידה.
- הפקודה `make valgrind` – בדיקת זליגת זיכרון באמצעות valgrind.
- הפקודה `make clean` - מוחקת את כל הקבצים הלא רלוונטיים לאחר ההרצה.


בהצלחה!


</div>
