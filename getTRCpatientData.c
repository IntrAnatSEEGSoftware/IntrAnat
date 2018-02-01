#include <stdio.h> 
#include <string.h>
/**
 *  ANONYMIZATION Of patient data in Micromed TRC files
 *
 * (c) Inserm U836 2014 - Manik Bhattacharjee
 *
 * License GNU GPL v3
*/

/**
  From TRC format description "System 98 - Type "4" EEG FILE STRUCTURE DESCRIPTION (Micromed)
  Offset of patient data, offset of acquisition date/time
*/
#define PATIENT_DATA_OFFSET 64
#define PATIENT_DATA_NAME_SIZE 22
#define PATIENT_DATA_FIRSTNAME_SIZE 20
#define PATIENT_BIRTH_DATE_SIZE 3
#define DATE_OFFSET 128
#define TIME_OFFSET 131
#define SPACE_CHAR 32 //Ascii code

int main(int argc, char **argv){
	int length, remains, i;
	char name[PATIENT_DATA_NAME_SIZE +1 ] = "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0";
	char firstname[PATIENT_DATA_FIRSTNAME_SIZE +1 ] = "\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0";
	unsigned char month, day, year, hours, minutes, seconds;
	if (argc != 2){
		printf("Usage : getTRCpatientData file.trc\n");
		return -1;
	}
	// Opening file read binary mode
	FILE *const fp = fopen(argv[1], "rb"); 
	if (!fp){
		printf("Cannot open TRC file\n");
		return -1;
	}
  if (fseek(fp, PATIENT_DATA_OFFSET, SEEK_SET) != 0){
		printf("Cannot go to Patient data section in TRC file !\n");
		return -1;
	}
	// Patient Last name
	fread (name, 1, PATIENT_DATA_NAME_SIZE, fp);
	
	//Patient First name
	fread (name, 1, PATIENT_DATA_FIRSTNAME_SIZE, fp);
	
	//Patient birth date
	month = fgetc(fp);
	day = fgetc(fp);
	year = fgetc(fp);
	printf("Subject : %s %s, born %d/%d/%d\n", firstname, name, day, month, year/100 > 0?2000+year%100:1900+year);
	
	if (fseek(fp, DATE_OFFSET, SEEK_SET) != 0){
		printf("Cannot go to Date data section in TRC file !\n");
		return -1;
	}
	day = fgetc(fp);
	month = fgetc(fp);
	year = fgetc(fp);
	
	if (fseek(fp, TIME_OFFSET, SEEK_SET) != 0){
		printf("Cannot go to Time data section in TRC file !\n");
		return -1;
	}
	hours = fgetc(fp);
	minutes = fgetc(fp);
	seconds = fgetc(fp);
	
	printf("Acquisition : %d/%d/%d %d:%d:%d\n", day, month, year/100 > 0?2000+year%100:1900+year, hours, minutes, seconds);
	printf("GRE_%d_%c%c%c%c\n", year/100 > 0?2000+year%100:1900+year, name[0], name[1], name[2], firstname[0]);
	
	fclose(fp);
}
