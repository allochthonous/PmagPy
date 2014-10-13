#!/usr/bin/env python

import wx
import controlled_vocabularies as vocab

# this module will provide all the functionality for the drop-down controlled vocabulary menus
# ideally, these classes should take in a grid item and then provide the functionality


class Menus():

    def __init__(self, data_type, check, grid, belongs_to, headers=['site_class', 'site_lithology', 'site_type', 'site_definition']):
        """take: data_type (string), check (top level class object for ErMagic steps 1-6), grid (grid object), belongs_to (options for data object to belong to, i.e. locations for the site Menus)"""
        self.data_type = data_type
        self.check = check # check is top level class object for entire ErMagic steps 1-6
        self.grid = grid
        self.window = grid.Parent # parent window in which grid resides
        self.belongs_to = belongs_to
        self.headers = headers
        self.selected_col = None
        self.InitUI()

    def InitUI(self):
        belongs_to = self.belongs_to
        if self.data_type == 'specimen':
            choices = {2: (belongs_to, False)}
        if self.data_type == 'sample' or self.data_type == 'site':
            choices = {2: (belongs_to, False), 3: (vocab.site_class, False), 4: (vocab.site_lithology, True), 5: (vocab.site_type, False)}
        if self.data_type == 'site':
            choices[6] = (vocab.site_definition, False)
        if self.data_type == 'location':
            choices = {1: (vocab.location_type, False)}
        if self.data_type == 'age':
            choices = {3: (vocab.geochronology_method_codes, False), 5: (vocab.age_units, False)}
        self.window.Bind(wx.grid.EVT_GRID_SELECT_CELL, lambda event: self.on_left_click(event, self.grid, choices), self.grid) 
        self.window.Bind(wx.grid.EVT_GRID_LABEL_LEFT_CLICK, self.on_label_click, self.grid)


    def on_label_click(self, event):
        col = event.GetCol()
        if col == -1:
            return 0
        if col == 2 and self.data_type == 'age':
            return 0
        if (col not in (-1, 0, 1)) or (col == 1 and self.data_type == 'location'):
        # if a new column was chosen without de-selecting the previous column, deselect the old selected_col
            if self.selected_col != None and self.selected_col != col: 
                col_label_value = self.grid.GetColLabelValue(self.selected_col)
                self.grid.SetColLabelValue(self.selected_col, col_label_value[:-10])
                for row in range(self.grid.GetNumberRows()):
                    self.grid.SetCellBackgroundColour(row, self.selected_col, 'white')
                self.grid.ForceRefresh()
            # deselect col if user is clicking on it a second time
            if col == self.selected_col:  
                col_label_value = self.grid.GetColLabelValue(col)
                self.grid.SetColLabelValue(col, col_label_value[:-10])
                for row in range(self.grid.GetNumberRows()):
                    self.grid.SetCellBackgroundColour(row, col, 'white')
                self.grid.ForceRefresh()
                self.selected_col = None
            else:
                self.selected_col = col
                col_label_value = self.grid.GetColLabelValue(col)
                self.grid.SetColLabelValue(col, col_label_value + " \nEDIT ALL")
                for row in range(self.grid.GetNumberRows()):
                    self.grid.SetCellBackgroundColour(row, col, 'light blue')
                self.grid.ForceRefresh()
        has_dropdown = ((col == 2 and self.data_type is 'specimen') or (col in range(2, 6) and self.data_type in ['site', 'sample']) or (col in (3, 5) and self.data_type == 'age') or (col == 6 and self.data_type == 'site') or (col == 1 and self.data_type == 'location'))

        # if the column has no drop-down list, allow user to edit all cells in the column through text entry
        if not has_dropdown and col not in (0, 1):  
            if self.selected_col == col:
                self.check.changes = True
                default_value = self.grid.GetCellValue(0, col)
                #data = wx.GetTextFromUser("Enter value for all cells in the column\nNote: this will overwrite any existing cell values", "Edit All", default_value)
                
                data = None
                dialog = wx.TextEntryDialog(None, "Enter value for all cells in the column\nNote: this will overwrite any existing cell values", "Edit All", default_value, style=wx.OK|wx.CANCEL)
                if dialog.ShowModal() == wx.ID_OK: 
                    data = dialog.GetValue() 
                    for row in range(self.grid.GetNumberRows()):
                        self.grid.SetCellValue(row, col, str(data))
                dialog.Destroy()
                # then deselect column
                col_label_value = self.grid.GetColLabelValue(col)
                self.grid.SetColLabelValue(col, col_label_value[:-10])
                for row in range(self.grid.GetNumberRows()):
                    self.grid.SetCellBackgroundColour(row, col, 'white')
                self.grid.ForceRefresh()
                self.selected_col = None


            
    def clean_up(self, grid):
        if self.selected_col:
            col_label_value = self.grid.GetColLabelValue(self.selected_col)
            self.grid.SetColLabelValue(self.selected_col, col_label_value[:-10])
            for row in range(self.grid.GetNumberRows()):
                self.grid.SetCellBackgroundColour(row, self.selected_col, 'white')
        self.grid.ForceRefresh()


    def on_left_click(self, event, grid, choices):
        """creates popup menu when user clicks on the column
        if that column is in the list of choices that get a drop-down menu.
        allows user to edit the column, but only from available values"""
        col = event.GetCol()
        if col in choices.keys(): # column should have a pop-up menu
            row = event.GetRow()
            menu = wx.Menu()
            two_tiered = choices[col][1]
            choices = choices[col][0]
            if not two_tiered: # menu is one tiered
                # insert CLEAR option as first choice
                if 'CLEAR cell of all values' not in choices:
                    choices.insert(0, 'CLEAR cell of all values')
                for choice in choices:
                    if not choice: choice = " " # prevents error if choice is an empty string
                    menuitem = menu.Append(wx.ID_ANY, choice)
                    self.window.Bind(wx.EVT_MENU, lambda event: self.on_select_menuitem(event, grid, row, col), menuitem)
                self.window.PopupMenu(menu)
                menu.Destroy()
            else: # menu is two_tiered
                # insert CLEAR option as first menuitem
                clear = menu.Append(-1, 'CLEAR cell of all values')
                self.window.Bind(wx.EVT_MENU, lambda event: self.on_select_menuitem(event, grid, row, col), clear)
                for choice in sorted(choices.items()):
                    submenu = wx.Menu()
                    for item in choice[1]:
                        menuitem = submenu.Append(-1, item)
                        self.window.Bind(wx.EVT_MENU, lambda event: self.on_select_menuitem(event, grid, row, col), menuitem)
                    menu.AppendMenu(-1, choice[0], submenu)
                self.window.PopupMenu(menu)
                menu.Destroy()

    def update_drop_down_menu(self, grid, choices):
        self.window.Bind(wx.grid.EVT_GRID_SELECT_CELL, lambda event: self.on_left_click(event, grid, choices), grid) 


    def on_select_menuitem(self, event, grid, row, col):
        """
        sets value of selected cell to value selected from menu
        """
        self.check.changes = True # if user selects a menuitem, that is an edit
        item_id =  event.GetId()
        item = event.EventObject.FindItemById(item_id)
        label = item.Label
        cell_value = grid.GetCellValue(row, col)
        if str(label) == "CLEAR cell of all values":
            label = ""
        elif (col in range(3, 6) and self.data_type in ['site', 'sample']) or (col == 3 and self.data_type == 'age'):
        #if True:
            if not label.lower() in cell_value.lower():
                label += (":" + cell_value).rstrip(':')
            else:
                label = cell_value
        if self.selected_col and self.selected_col == col:
            for row in range(self.grid.GetNumberRows()):
                grid.SetCellValue(row, col, label)
                #self.selected_col = None
        else:
            grid.SetCellValue(row, col, label)


