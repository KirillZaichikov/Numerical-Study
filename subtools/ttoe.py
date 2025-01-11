import openpyxl
from openpyxl.styles import Font

def dat_to_xlsx(dat_filepath, xlsx_filepath):
    try:
        with open(dat_filepath, 'r') as dat_file:
            data = []
            for line in dat_file:
                try:
                    x, y, *rest = map(float, line.strip().split())
                    data.append([x, y])
                except ValueError:
                    print(f"Ошибка при обработке строки: {line.strip()}. Проверьте формат данных.")

        workbook = openpyxl.Workbook()
        sheet = workbook.active
        sheet.append(["X", "Y"])

        for row in data:
            sheet.append(row)

        workbook.save(xlsx_filepath)
        print(f"Данные успешно записаны в файл {xlsx_filepath}")

    except FileNotFoundError:
        print(f"Файл {dat_filepath} не найден.")
    except Exception as e:
        print(f"Произошла ошибка: {e}")


if __name__ == '__main__':
    dat_file_path = "TRNN.dat"
    xlsx_file_path = "TRNN.xlsx"
    dat_to_xlsx(dat_file_path, xlsx_file_path)