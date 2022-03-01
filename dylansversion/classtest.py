class Employee:
    def __init__(self, first, last):
        self.first = first
        self.last = last
    
    def fullname(self):
        return f'{self.first} {self.last}'


emp_1 = Employee('Dylan', 'Jones')
emp_2 = Employee('Cooper', 'Ritter')



print(emp_2.fullname())

