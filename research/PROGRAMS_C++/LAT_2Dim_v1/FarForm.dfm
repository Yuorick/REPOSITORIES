object Form2: TForm2
  Left = 0
  Top = 0
  Caption = #1060#1040#1056
  ClientHeight = 489
  ClientWidth = 784
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -19
  Font.Name = 'Tahoma'
  Font.Style = [fsBold]
  OldCreateOrder = False
  Position = poOwnerFormCenter
  OnShow = FormShow
  PixelsPerInch = 96
  TextHeight = 23
  object Label1: TLabel
    Left = 136
    Top = 36
    Width = 42
    Height = 23
    Caption = #1060#1040#1056
  end
  object Button1: TButton
    Left = 112
    Top = 415
    Width = 75
    Height = 25
    Caption = 'OK'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -16
    Font.Name = 'MS Sans Serif'
    Font.Style = []
    ParentFont = False
    TabOrder = 0
    OnClick = Button1Click
  end
  object Panel1: TPanel
    Left = 8
    Top = 80
    Width = 345
    Height = 281
    BevelOuter = bvNone
    TabOrder = 1
    object Image1: TImage
      Left = 0
      Top = 0
      Width = 345
      Height = 281
      Align = alClient
      OnMouseDown = Image1MouseDown
      ExplicitLeft = 40
      ExplicitTop = 104
      ExplicitWidth = 281
      ExplicitHeight = 100
    end
  end
end
