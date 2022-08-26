object TkachenkoForm: TTkachenkoForm
  Left = 0
  Top = 0
  Caption = #1056#1072#1089#1089#1095#1077#1090' '#1076#1083#1103' '#1058#1082#1072#1095#1077#1085#1082#1086
  ClientHeight = 237
  ClientWidth = 792
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  Position = poMainFormCenter
  PixelsPerInch = 96
  TextHeight = 13
  object Label1: TLabel
    Left = 8
    Top = 8
    Width = 179
    Height = 18
    Caption = #1044#1080#1088#1077#1082#1090#1086#1088#1080#1103' '#1088#1077#1079#1091#1083#1100#1090#1072#1090#1086#1074
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -16
    Font.Name = #1060#1082#1096#1092#1076
    Font.Style = []
    ParentFont = False
  end
  object Label2: TLabel
    Left = 16
    Top = 144
    Width = 227
    Height = 18
    Caption = #1042#1086#1079#1084#1086#1078#1085#1099#1077' '#1084#1086#1097#1085#1086#1089#1090#1080' '#1085#1072' '#1082#1072#1085#1072#1083
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -16
    Font.Name = 'Arial'
    Font.Style = []
    ParentFont = False
  end
  object Label3: TLabel
    Left = 16
    Top = 80
    Width = 171
    Height = 18
    Caption = #1042#1086#1079#1084#1086#1078#1085#1099#1077' '#1090#1080#1087#1099' '#1094#1077#1083#1077#1081
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -16
    Font.Name = 'Arial'
    Font.Style = []
    ParentFont = False
  end
  object Label4: TLabel
    Left = 368
    Top = 144
    Width = 157
    Height = 18
    Caption = #1055#1077#1088#1077#1095#1077#1085#1100' '#1074#1099#1089#1086#1090' '#1094#1077#1083#1080
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -16
    Font.Name = 'Arial'
    Font.Style = []
    ParentFont = False
  end
  object Label5: TLabel
    Left = 368
    Top = 80
    Width = 200
    Height = 18
    Caption = #1042#1086#1079#1084#1086#1078#1085#1072#1103' '#1075#1077#1086#1084#1077#1090#1088#1080#1103' '#1060#1040#1056
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -16
    Font.Name = 'Arial'
    Font.Style = []
    ParentFont = False
  end
  object Edit1: TEdit
    Left = 8
    Top = 32
    Width = 737
    Height = 27
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -16
    Font.Name = 'Tahoma'
    Font.Style = []
    ParentFont = False
    TabOrder = 0
    Text = 'Edit1'
  end
  object Button1: TButton
    Left = 751
    Top = 32
    Width = 33
    Height = 27
    Caption = '...'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -21
    Font.Name = 'Arial'
    Font.Style = [fsBold]
    ParentFont = False
    TabOrder = 1
    OnClick = Button1Click
  end
  object CB3: TComboBox
    Left = 255
    Top = 141
    Width = 50
    Height = 26
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -16
    Font.Name = 'Arial'
    Font.Style = []
    ItemIndex = 0
    ParentFont = False
    TabOrder = 2
    Text = '1'
    Items.Strings = (
      '1'
      '5'
      '10')
  end
  object Button2: TButton
    Left = 8
    Top = 192
    Width = 144
    Height = 25
    Caption = #1056#1040#1057#1057#1063#1048#1058#1040#1058#1068
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -16
    Font.Name = 'Arial'
    Font.Style = []
    ParentFont = False
    TabOrder = 3
    OnClick = Button2Click
  end
  object CB1: TComboBox
    Left = 193
    Top = 77
    Width = 120
    Height = 26
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -16
    Font.Name = 'Arial'
    Font.Style = []
    ItemIndex = 0
    ParentFont = False
    TabOrder = 4
    Text = #1043#1072#1088#1087#1091#1085'300'
    Items.Strings = (
      #1043#1072#1088#1087#1091#1085'300'
      #1043#1072#1088#1087#1091#1085'700'
      #1057#1072#1084#1086#1083#1105#1090)
  end
  object CB4: TComboBox
    Left = 545
    Top = 141
    Width = 64
    Height = 26
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -16
    Font.Name = 'Arial'
    Font.Style = []
    ItemIndex = 0
    ParentFont = False
    TabOrder = 5
    Text = '5'
    Items.Strings = (
      '5'
      '10'
      '20'
      '50'
      '100'
      '300'
      '500')
  end
  object CB2: TComboBox
    Left = 583
    Top = 77
    Width = 59
    Height = 26
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -16
    Font.Name = 'Arial'
    Font.Style = []
    ParentFont = False
    TabOrder = 6
    Text = '8x8'
    Items.Strings = (
      '8x8'
      '4x8'
      '4x12'
      '4x16'
      '4x20'
      '6x12')
  end
  object OpenDlg: TOpenDialog
    Left = 656
    Top = 184
  end
end
