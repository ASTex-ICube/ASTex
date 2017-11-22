/*******************************************************************************
* ASTex:                                                                       *
* Copyright (C) IGG Group, ICube, University of Strasbourg, France             *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: https://astex-icube.github.io                                      *
* Contact information: astex@icube.unistra.fr                                  *
*                                                                              *
*******************************************************************************/




#ifndef IMAGEVIEWER_H
#define IMAGEVIEWER_H

#include <QApplication>
#include <QMainWindow>
#include <QImage>
#include <QtWidgets>
#include <QFileDialog>

#include <cstring>
#include <iostream>
#include <ASTex/internal.h>

QT_BEGIN_NAMESPACE
class QAction;
class QLabel;
class QMenu;
class QScrollArea;
class QScrollBar;
QT_END_NAMESPACE

#include<ASTex/image_common.h>

// forward (must be define whehe .h included
void app_mouse_clicked(int button, int x, int y,int id);
void app_key_pressed(int code, char key, int);

class ImageViewer : public QMainWindow
{
	Q_OBJECT

public:
	/**
	 * @brief Constructor
	 * @param name windows name
	 * @param app appliction ptr (for out with ESC)
	 * @param id id of win (for use of app_mouse_clicked & app_key_pressed)
	 */
	ImageViewer(const std::string& name, QApplication* app=nullptr, int id=0);
	
	template <typename IMG>
	inline auto update(const IMG& img, int zoom=1) -> typename std::enable_if<IMG::NB_CHANNELS==1>::type
	{
		using SCALAR = typename IMG::DataType;

		const SCALAR* ptr = img.getDataPtr();
		int w = img.width();
		int h = img.height();

		auto conv = [&] (const SCALAR x) -> uint8_t
		{
			if (std::is_floating_point<SCALAR>::value)
				return ASTex::unnormalized<uint8_t>(x);
			return uint8_t(x);
		};

		QImage im(w,h,QImage::Format_RGB888);
		for(int i=0;i<h;++i)
		{
			uint8_t *optr = im.scanLine(i);
			for(int j=0;j<w;++j)
			{
				*optr++ = conv(*ptr);
				*optr++ = conv(*ptr);
				*optr++ = conv(*ptr++);
			}
		}

		if (zoom != 1)
		{
			QImage imz = im.scaledToWidth(w*zoom);
			im.swap(imz);
		}

		m_imageLabel->setPixmap(QPixmap::fromImage(im));
		if ((w<=1024) && (h<=1024))
			m_scrollArea->setMinimumSize(im.width()+2, im.height()+2);

		m_scrollArea->setVisible(true);
		m_imageLabel->adjustSize();
	}

	template <typename IMG>
	inline auto update(const IMG& img, int zoom=1) -> typename std::enable_if<IMG::NB_CHANNELS==3>::type
	{
		using SCALAR = typename IMG::DataType;

		const SCALAR* ptr = img.getDataPtr();
		int w = img.width();
		int h = img.height();

		auto conv = [&] (const SCALAR x) -> uint8_t
		{
			if (std::is_floating_point<SCALAR>::value)
				return ASTex::unnormalized<uint8_t>(x);
			return uint8_t(x);
		};

		QImage im(w,h,QImage::Format_RGB888);
		for(int i=0;i<h;++i)
		{
			uint8_t *optr = im.scanLine(i);
			for(int j=0;j<w;++j)
			{
				*optr++ = conv(*ptr++);
				*optr++ = conv(*ptr++);
				*optr++ = conv(*ptr++);
			}
		}
		if (zoom != 1)
		{
			QImage imz = im.scaledToWidth(w*zoom);
			im.swap(imz);
		}

		m_imageLabel->setPixmap(QPixmap::fromImage(im));
		if ((w<=1024) && (h<=1024))
			m_scrollArea->setMinimumSize(im.width()+2, im.height()+2);

		m_scrollArea->setVisible(true);
		m_imageLabel->adjustSize();
	}

	template <typename IMG>
	inline auto update(const IMG& img, int zoom=1) -> typename std::enable_if<IMG::NB_CHANNELS==4>::type
	{
		using SCALAR = typename IMG::DataType;

		const SCALAR* ptr = img.getDataPtr();
		int w = img.width();
		int h = img.height();

		auto conv = [&] (const SCALAR x) -> uint8_t
		{
			if (std::is_floating_point<SCALAR>::value)
				return ASTex::unnormalized<uint8_t>(x);
			return uint8_t(x);
		};

		QImage im(w,h,QImage::Format_RGB888);
		for(int i=0;i<h;++i)
		{
			uint8_t *optr = im.scanLine(i);
			for(int j=0;j<w;++j)
			{
				*optr++ = conv(*ptr++);
				*optr++ = conv(*ptr++);
				*optr++ = conv(*ptr++);
				*optr++ = conv(*ptr++);
			}
		}
		if (zoom != 1)
		{
			QImage imz = im.scaledToWidth(w*zoom);
			im.swap(imz);
		}

		m_imageLabel->setPixmap(QPixmap::fromImage(im));
		if ((w<=1024) && (h<=1024))
			m_scrollArea->setMinimumSize(im.width()+2, im.height()+2);

		m_scrollArea->setVisible(true);
		m_imageLabel->adjustSize();
	}

	/**
	 * get horizontal scroll bar for easy connect/syncro
	 */
	inline const QScrollBar* hsb() const { return m_scrollArea->horizontalScrollBar(); }
	
	/**
	 * get vertical scroll bar for easy connect/syncro
	 */
	inline const QScrollBar* vsb() const { return m_scrollArea->verticalScrollBar(); }

private:
	QLabel* m_imageLabel;
	QScrollArea* m_scrollArea;
	int m_id;

	inline void mousePressEvent(QMouseEvent* event)
	{
		app_mouse_clicked(event->button(),
						  event->x()-m_scrollArea->x()+m_scrollArea->horizontalScrollBar()->value(),
						  event->y()-m_scrollArea->y()+m_scrollArea->verticalScrollBar()->value(),
						  m_id);
	}

	inline void keyPressEvent(QKeyEvent* event)
	{
		app_key_pressed(event->key(),event->text().toStdString()[0],m_id);
	}


	inline void resizeEvent(QResizeEvent* event)
	{
		QMainWindow::resizeEvent(event);
		emit(win_resized(event->size()));
	}

signals:
   void win_resized(QSize sz);

public slots:
	inline void win_resize(QSize sz) { resize(sz); }
};


inline ImageViewer::ImageViewer(const std::string& name, QApplication* app, int id):
	m_imageLabel(new QLabel), m_scrollArea(new QScrollArea),m_id(id)
{
	m_imageLabel->setBackgroundRole(QPalette::Base);
	m_imageLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
	m_imageLabel->setScaledContents(true);
	m_scrollArea->setBackgroundRole(QPalette::Dark);
	m_scrollArea->setWidget(m_imageLabel);
	m_scrollArea->setVisible(false);
	setCentralWidget(m_scrollArea);

	if (app != nullptr)
	{
		QMenu* fileMenu = menuBar()->addMenu(tr("&File"));
		QAction *exitAct = fileMenu->addAction(tr("E&xit"), app, SLOT(closeAllWindows()));
		exitAct->setShortcut(tr("Esc"));
	}

	setWindowTitle(QString(name.c_str()));
}



//template <typename IMG>
//std::unique_ptr<ImageViewer> image_viewer(const IMG& img, const std::string& name="", QApplication* app=nullptr, int id=0)
//{
//	std::unique_ptr<ImageViewer> view(new ImageViewer(name,app,id));
//	view->update(img);
//	view->show();
//}

template <typename IMG>
inline std::unique_ptr<ImageViewer> image_viewer(const IMG& img, const std::string& name="", QApplication* app=nullptr, int id=0)
{
	ImageViewer* view = new ImageViewer(name,app,id);
	view->update(img);
	view->show();
	return std::unique_ptr<ImageViewer>(view);
}

#endif	
